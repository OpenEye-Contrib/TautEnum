//
// file taut_enum.cc
// David Cosgrove
// AstraZeneca
// 31st August 2011
//
// This is a standalone program that enumerates all tautomers of the molecules
// input. It does this in two stages - first it standardises them, then
// it enumerates from the standardised set.  To do this, it can be passed
// a SMIRKS file for each or either step, but by default it has a set of
// canned SMIRKS definitions to use, built in in the class TautStand and
// TautEnum.  There are 2 levels of enumeration. One corresponds at this
// point to Pete Kenny's Leatherface enumerations, the other is much more
// ambitious but may well produce tautomers that aren't at all likely to
// be seen in real life.
//

#include "TautEnum.H"
#include "TautEnumCallableBase.H"
#include "TautEnumCallableSerial.H"
#include "TautEnumCallableThreaded.H"
#include "TautEnumSettings.H"
#include "TautStand.H"
#include "FileExceptions.H"

#include <iostream>
#include <list>

#include <oechem.h>

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/thread.hpp>

using namespace std;
using namespace OEChem;
using namespace OESystem;

namespace DACLIB {
void apply_daylight_aromatic_model( OEMolBase &mol );
}

extern string BUILD_TIME; // in build_time.cc

typedef boost::shared_ptr<OEMolBase> pOEMol;

// ****************************************************************************
// fix SMILES so AtomMap indices don't appear. A spell from James H. at Openeye
void fix_output_smiles_format( oemolstreambase &oms ) {

  // both of these formats need to write Canonical SMILES strings which is not
  // default behaviour
  if( oms.GetFormat() == OEFormat::SMI ) {
    unsigned int oflavor = OEOFlavor::SMI::Default;
    oflavor = oflavor ^ OEOFlavor::SMI::AtomMaps;
    oflavor = oflavor | OEOFlavor::SMI::Canonical;
    oms.SetFlavor( OEFormat::SMI , oflavor );
  } else if( oms.GetFormat() == OEFormat::ISM ) {
    unsigned int oflavor = OEOFlavor::ISM::Default^OEOFlavor::SMI::AtomMaps;
    oflavor = oflavor | OEOFlavor::SMI::Canonical | OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo ;
    oms.SetFlavor( OEFormat::ISM, oflavor );
  }

}

// ****************************************************************************
void serial_run( const TautEnumSettings &tes ) {

  cerr << "Serial run." << endl;

  oemolistream ims;
  if( !ims.open( tes.input_mol_file() ) ) {
    cerr << "Failed to open " << tes.input_mol_file() << " for reading." << endl;
    exit( 1 );
  }

  oemolostream oms;
  if( !oms.open( tes.output_mol_file() ) ) {
    cerr << "Failed to open " << tes.output_mol_file() << " for writing." << endl;
    exit( 1 );
  }

  fix_output_smiles_format( oms );

  TautEnumCallableSerial tc( &ims , &oms , tes );

  tc();

}

// ****************************************************************************
void threaded_run( const TautEnumSettings &tes ) {

  cerr << endl
       << "Currently, for reasons which I don't understand, the threaded mode" << endl
       << "takes longer and is less reliable than non-threaded.  It is wont to" << endl
       << "hang consuming CPU but without writing further output, for example." << endl
       << "Its use is not currently recommended." << endl << endl;

  // set up the threaded input and output streams
  boost::thread_group tg;
  oemolithread ims;
  ims.open( tes.input_mol_file() );

  oemolothread oms;
  oms.open( tes.output_mol_file() );
  fix_output_smiles_format( oms );

  // In OEToolkits 1.7.6, OEPerceiveChiral, which is used in TautEnum, gives a memory error
  // using the default memory pool system.  Either of these two fixes it, at the expense of
  // speed.  The 2nd one seems to be a bit faster.
  // OESystem::OESetMemPoolMode( OESystem::OEMemPoolMode::System );
  OESystem::OESetMemPoolMode(OESystem::OEMemPoolMode::Mutexed|OESystem::OEMemPoolMode::UnboundedCache);

  TautEnumCallableThreaded tct( &ims , &oms , tes );

  // create the threads
  list<TautEnumCallableThreaded> callables;
  int nt = 1;
  if( tes.num_threads() <= 0 ) {
    nt = boost::thread::hardware_concurrency() + tes.num_threads();
  } else {
    nt = tes.num_threads();
  }

  cerr << "Number of threads to use : " << nt << endl;

  for( int i = 0 ; i < nt ; i++ ) {
    callables.push_back( tct );
    tg.create_thread( boost::ref( callables.back() ) ); // careful not to copy callable again
  }

  tg.join_all(); // wait for them all to finish
  boost::this_thread::sleep( boost::posix_time::seconds( 1 ) ); // give worker threads time to deallocate

}

// ****************************************************************************
int main( int argc , char **argv ) {

  cerr << endl << "taut_enum, built " << BUILD_TIME << " using OEToolkits version "
       << OEChem::OEChemGetRelease() << " (" << OEChem::OEChemGetVersion() << ")." << endl;

  TautEnumSettings tes( argc , argv );

  if( !tes ) {
    tes.print_error( cout );
    tes.print_error( cerr );
    tes.print_usage( cout );
    exit( 1 );
  }

  // Suppress irritating warnings from OELibraryGen (and everything else, of course, but
  // OELibraryGen gives a lot of very irritating stuff)
  OESystem::OEThrow.SetLevel( OESystem::OEErrorLevel::Error );

  if( tes.do_threaded() ) {
    threaded_run( tes );
  } else {
    serial_run( tes );
  }

}
