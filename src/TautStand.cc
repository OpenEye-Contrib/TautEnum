//
// file TautStand.cc
// David Cosgrove
// AstraZeneca
// 2nd August 2011
//

#include "TautStand.H"

#include <iostream>

#include <oechem.h>
#include <oeplatform.h>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

typedef boost::shared_ptr<OEChem::OELibraryGen> pOELibGen;
typedef boost::shared_ptr<OEChem::OEMolBase> pOEMolBase;

// ****************************************************************************
// in smirks_helper_fns.cc
namespace DACLIB {
string create_cansmi( const OEMolBase &in_mol );
void read_vbs_from_file( const string &filename ,
                         vector<pair<string,string> > &vbs );
void read_vbs_from_string( const string &vbs_string , vector<pair<string,string> > &vbs );
void read_smirks_from_file( const string &filename ,
                            vector<pair<string,string> > &smks );
void read_smirks_from_string( const string &smks_string , vector<pair<string,string> > &smks );
void expand_vector_bindings( const vector<pair<string,string> > &in_smirks ,
                             vector<pair<string,string> > &vbs ,
                             vector<string> &exp_smirks );
void create_libgens( const vector<string> &exp_smirks ,
                     const vector<pair<string,string> > &in_smirks ,
                     vector<pOELibGen> &lib_gens );
}

// ****************************************************************************
TautStand::TautStand( const string &smirks_string , const string &vb_string ) {

  DACLIB::read_vbs_from_string( vb_string , vbs_ );
  DACLIB::read_smirks_from_string( smirks_string , smirks_ );

  DACLIB::expand_vector_bindings( smirks_ , vbs_ , exp_smirks_ );

}

// ****************************************************************************
TautStand::TautStand( const string &smirks_file , const string &vb_file ,
                      bool dummy __attribute__((unused)) ) :
  smirks_file_( smirks_file ) , vb_file_( vb_file ) {

#ifdef NOTYET
  cout << "loading standardisation smirks from " << smirks_file
       << " with vbs from " << vb_file << endl;
#endif

  DACLIB::read_smirks_from_file( smirks_file_ , smirks_ );
  DACLIB::read_vbs_from_file( vb_file , vbs_ );
  DACLIB::expand_vector_bindings( smirks_ , vbs_ , exp_smirks_ );

}

// ****************************************************************************
// copy c'tor, needed for threading.
TautStand::TautStand( const TautStand &rhs ) {

  smirks_file_ = rhs.smirks_file_;
  vb_file_ = rhs.vb_file_;
  smirks_ = rhs.smirks_;
  vbs_ = rhs.vbs_;
  exp_smirks_ = rhs.exp_smirks_;

  // lib_gens_ is filled by standardise as required, so not copying it here. A deep
  // copy would have been required otherwise, I mention for future reference.

}

// ****************************************************************************
OEMolBase *TautStand::standardise( OEMolBase &in_mol , bool verbose ,
                                   bool add_smirks_to_name , bool strip_salts ) {

#ifdef NOTYET
  cout << "Standardising " << DACLIB::create_cansmi( in_mol )
       << " strip_salts : " << strip_salts << endl;
#endif

  // make the libgen objects up front if not already done
  if( lib_gens_.empty() ) {
    DACLIB::create_libgens( exp_smirks_ , smirks_ , lib_gens_ );
  }

  pOEMolBase prod_mol( OENewMolBase( in_mol , OEMolBaseType::OEDefault ) );
  // keep track of all intermediate SMILES strings in case we go round in an infinite loop.
  // most likely that will be a tautomer flipping backwards and forwards, but in principle it
  // could be a loop of more tautomers.  The SMIRKS aren't supposed to create such loops
  // but you should never underestimate the ability of a chemist to screw you over.
  // If this happens, return the last one found. It's in a standard form, after all, so should
  // be fine for further use.
  set<string> all_smis;
  all_smis.insert( DACLIB::create_cansmi( in_mol ) );

  while( true ) {
    unsigned int smis_size = all_smis.size();
    int smirks_num = 0;
    BOOST_FOREACH( pOELibGen libgen , lib_gens_ ) {
#ifdef NOTYET
      cout << "Next SMIRKS " << smirks_[smirks_num].first << " : " << smirks_[smirks_num].second << endl;
#endif
      libgen->SetAssignMapIdx( false ); // don't want them showing for this
      while( 1 ) {
        // SetStartingMaterial returns the number of matches of the SMIRKS in
        // the molecule, so don't do anything if it returns 0
        if( libgen->SetStartingMaterial( *prod_mol , 0 , false ) ) {
          OEIter<OEMolBase> prod = libgen->GetProducts();
          prod_mol.reset( OENewMolBase( *prod , OEMolBaseType::OEDefault ) );
          if( strip_salts ) {
            OETheFunctionFormerlyKnownAsStripSalts( *prod_mol );
          }
          OEFindRingAtomsAndBonds( *prod_mol );
          OEAssignAromaticFlags( *prod_mol );
          OEPerceiveChiral( *prod_mol );
          string this_smi = DACLIB::create_cansmi( *prod_mol );
          if( !all_smis.insert( this_smi ).second ) {
            cerr << "Problem with TautStand : " << in_mol.GetTitle()
                 << " creates an infinite loop of tautomers." << endl;
            break;
          }
          if( verbose ) {
            cout << "New product in tautomer standardiser : " << DACLIB::create_cansmi( *prod_mol ) << endl
                 << "Made from " << smirks_[smirks_num].first << " : " << smirks_[smirks_num].second << endl;
          }
          if( add_smirks_to_name ) {
            string curr_name = prod_mol->GetTitle();
            curr_name += string( " " ) + smirks_[smirks_num].first;
            prod_mol->SetTitle( curr_name );
          }
        } else {
          break;
        }
      }
      ++smirks_num;
    }
    if( all_smis.size() == smis_size ) {
      break; // didn't add anything new
    }
  }
#ifdef NOTYET
  cout << "Final answer : " << DACLIB::create_cansmi( *prod_mol ) << endl;
#endif

  return OENewMolBase( *prod_mol , OEMolBaseType::OEDefault );

}
