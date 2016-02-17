//
// File TautEnumCallable.H
// David Cosgrove
// AstraZeneca
// 8th February 2012.
//
// This public class is used to do threaded tautomer enumeration
// using boost threads.  It wraps a load of stuff that used to be in
// the taut_enum.cc main program so that it can run enumerations
// concurrently.  The output order of the molecules is not guaranteed
// to be the same as the input order, although all tautomers of the same
// molecule will be together.
//
// This class relies very heavily on work by Mikko Vainio who kindly
// came over from Molndal and showed me how to do it.

#ifndef TAUTENUMCALLABLE_H
#define TAUTENUMCALLABLE_H

#include "TautEnumSettings.H"

#include <string>

namespace OEChem {

class oemolithread;
class oemolothread;
class oemolstreambase;

}

class TautStand;
class TautEnum;

// ****************************************************************************

class TautEnumCallableBase {

public :

  TautEnumCallableBase( const TautEnumSettings &settings ) :
    tes_( settings ) {}
  TautEnumCallableBase( const TautEnumCallableBase &rhs ) : tes_( rhs.tes_ ) {}

  virtual ~TautEnumCallableBase() {};

  virtual void operator()(); // the operator that boost::thread calls to do the work

protected :

  TautEnumSettings tes_;

  // make the TautStand and TautEnum objects, using the relevant data from tes_
  virtual void create_enumerator_objects( const std::string &stand_smirks_file ,
                                          const std::string &enum_smirks_file ,
                                          const std::string &vb_file ,
                                          const std::string &default_stand_smirks ,
                                          const std::string &default_enum_smirks ,
                                          const std::string &default_vbs ,
                                          TautStand *&taut_stand , TautEnum *&taut_enum );

  virtual bool read_next_molecule( OEChem::OEMolBase &mol ) = 0;
  virtual void write_molecule( OEChem::OEMolBase &mol ) = 0;
  virtual void output_molecules( std::vector<OEChem::OEMolBase *> &out_mols );

};

void sort_and_uniquify_molecules( std::vector<OEChem::OEMolBase *> &mols );

#endif // TAUTENUMCALLABLE_H
