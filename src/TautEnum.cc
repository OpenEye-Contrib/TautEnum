//
// file TautEnum.cc
// David Cosgrove
// AstraZeneca
// 3rd August 2011
//

#include "TautEnum.H"
#include "chrono.h"

#include <oechem.h>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
// #include <boost/thread.hpp>

using namespace boost;
using namespace std;
using namespace OEChem;
using namespace OESystem;

// ****************************************************************************
// in smirks_helper_fns.cc
namespace DACLIB {
void read_vbs_from_file( const string &filename , vector<pair<string,string> > &vbs );
void read_vbs_from_string( const string &vbs_string , vector<pair<string,string> > &vbs );
void read_smirks_from_file( const string &filename , vector<pair<string,string> > &smks );
void read_smirks_from_string( const string &smks_string , vector<pair<string,string> > &smks );
void expand_vector_bindings( const vector<pair<string,string> > &in_smirks ,
                             vector<pair<string,string> > &vbs ,
                             vector<string> &exp_smirks );
string create_cansmi( const OEMolBase &in_mol );
void create_libgens( const vector<string> &exp_smirks ,
                     const vector<pair<string,string> > &in_smirks ,
                     vector<pOELibGen> &lib_gens );
OESubSearch *create_oesubsearch( const string &smarts , bool reorder ); // in eponymous file
string extract_smarts_from_smirks( const string &smirks ); // in eponymous file
void radical_atoms( OEMolBase &mol , vector<OEAtomBase *> &rad_atoms ); // in eponymous file
}

// ****************************************************************************
// in canned_tautenum_routines.cc
OEMolBase *build_copy_of_mol( OEMolBase &mol );

// ****************************************************************************
// 1 and only 1 of original_enumeration or extended_enumeration must be true
TautEnum::TautEnum( const string &smirks_string , const string &vbs_string ,
                    unsigned int max_t ) : max_out_mols_( max_t ) {

#ifdef NOTYET
  cout << "Loading enumeration SMIRKS from string" << endl << smirks_string << endl;
#endif

  DACLIB::read_vbs_from_string( vbs_string , vbs_ );
  DACLIB::read_smirks_from_string( smirks_string , smirks_ );
  DACLIB::expand_vector_bindings( smirks_ , vbs_ , exp_smirks_ );

#ifdef NOTYET
  cout << "Number of expanded enumeration SMIRKS : " << exp_smirks_.size() << endl;
#endif

}

// ****************************************************************************
TautEnum::TautEnum( const string &smirks_file , const string &vb_file ,
                    bool dummy __attribute__((unused)) , unsigned int max_t ) :
  smirks_file_( smirks_file ) , vb_file_( vb_file ) , max_out_mols_( max_t ) {

#ifdef NOTYET
  cout << "loading enumeration smirks from " << smirks_file
       << " with vbs from " << vb_file << endl;
#endif

  DACLIB::read_smirks_from_file( smirks_file_ , smirks_);
  DACLIB::read_vbs_from_file( vb_file , vbs_ );
  DACLIB::expand_vector_bindings( smirks_ , vbs_ , exp_smirks_ );

#ifdef NOTYET
  cout << "Number of expanded enumeration SMIRKS : " << exp_smirks_.size() << endl;
#endif

}

// ****************************************************************************
// copy c'tor, needed for threading.
TautEnum::TautEnum( const TautEnum &rhs ) : max_out_mols_( rhs.max_out_mols_ ) {

  smirks_file_ = rhs.smirks_file_;
  vb_file_ = rhs.vb_file_;
  smirks_ = rhs.smirks_;
  vbs_ = rhs.vbs_;
  exp_smirks_ = rhs.exp_smirks_;

  // lib_gens_ is filled from exp_smirks_ as required, so not copying it here.  A deep
  // copy would have been required otherwise, I mention for future reference.

}

// ****************************************************************************
// throws an exception of type TooManyTautomers if max_tauts_ is exceeded
// The input tautomer is always returned first in the output vector.
vector<OEMolBase *> TautEnum::enumerate( OEMolBase &in_mol , bool verbose ,
                                         bool add_smirks_to_name ) {

#ifdef NOTYET
  cout << "Generating tautomers for " << in_mol.GetTitle() << " : "
       << DACLIB::create_cansmi( in_mol ) << endl;
#endif

  set<string> all_can_smis;
  all_can_smis.insert( DACLIB::create_cansmi( in_mol ) );

  vector<OEMolBase *> ret_mols;
  ret_mols.push_back( OENewMolBase( in_mol , OEMolBaseType::OEDefault ) );

  // make the libgen objects up front if not already done
  if( lib_gens_.empty() ) {
    DACLIB::create_libgens( exp_smirks_ , smirks_ , lib_gens_ );
  }

  vector<OEAtomBase *> input_rad_atoms;
  DACLIB::radical_atoms( in_mol , input_rad_atoms );

  size_t next_start = 0;
  while( true ) {
    // only do tautomers added in the last round. There should be no further products
    // possible from the results of rounds previous to that, as they will already be
    // in ret_mols.
    size_t start_size = ret_mols.size();
#ifdef NOTYET
    cout << "Next start, current set are : " << endl;
    BOOST_FOREACH( string smi , all_can_smis ) {
      cout << smi << endl;
    }
#endif

    for( size_t i = next_start , is = ret_mols.size() ; i < is ; ++i ) {

      int smirks_num = 0;
      BOOST_FOREACH( pOELibGen libgen , lib_gens_ ) {

#ifdef NOTYET
        cout << "NEXT SMIRKS : " << smirks_[smirks_num].first << " : " << smirks_[smirks_num].second
             << " : " << exp_smirks_[smirks_num] << endl << endl;
#endif

        libgen->SetAssignMapIdx( false ); // don't need map indices for this
        libgen->SetStartingMaterial( *ret_mols[i] , 0 , false );
        // this is a new function from 2013.Feb beta release that we're testing
        // at the moment.
        libgen->SetValidateKekule( false );

        OEIter<OEMolBase> prod = libgen->GetProducts();
        if( prod ) {
#ifdef NOTYET
          // this isn't really needed any more.  Run taut_enum with --verbose.
          cout << endl << "Prod for next libgen" << endl;
          cout << "SMIRKS : " << smirks_[smirks_num].first << " : " << smirks_[smirks_num].second << endl;
          cout << "Expanded SMIRKS : " << exp_smirks_[smirks_num] << endl;
          string inputsmi;
          OECreateCanSmiString( inputsmi , *ret_mols[i] );
          cout << "Input SMILES : " << inputsmi << endl;
#endif
          for( ; prod ; ++prod ) {
#ifdef NOTYET
            cout << "raw prod_mol : " << DACLIB::create_cansmi( *prod ) << endl;
#endif
            // Up to OEToolkits v 2012.Oct (v1.9.0) some molecules with extended
            // aromaticity got screwed up by some of the SMIRKS. e.g.
            // c1ccc2c(c1)c(=O)c3ccc4c(c3c2=O)[nH]c5ccc6c(=O)ccc(=O)c6c5[nH]4
            // when tackled with
            // SMIRKS : ENUM_AROM_9_4 : [H:8][n;H1;X3;!+:1]:[$CAR:2]:[$CAR:3]:[$CAR:4]:[$CAR:5]:[c:6]=[$REV_OS:7]>>[*:1]:[*:2]:[*:3]:[*:4]:[*:5]:[c:6]-[*:7][H:8]
            // gives, inter alia, c1ccc2c(c1)C(=O)c3ccc4c(c3C2=O)nc5ccc6c(c5n4)C(=O)[CH]C=C6O
            // where similar rings such as c1cc2c(c3c1[nH]c4c5c(cc(c4[nH]3))c(=O)c6ccccc6c5=O)c(=O)c7ccccc7c2=O
            // are ok.
            OEMolBase *prod_mol = OENewMolBase( *prod , OEMolBaseType::OEDefault );
            OEFindRingAtomsAndBonds( *prod_mol );
            OEAssignAromaticFlags( *prod_mol );
            OEPerceiveChiral( *prod_mol );
            vector<OEAtomBase *> prod_rad_atoms;
            DACLIB::radical_atoms( *prod_mol , prod_rad_atoms );
            if( prod_rad_atoms.size() > input_rad_atoms.size() ) {
              // we don't want products that have created free radicals. We'd rather the SMIRKS
              // toolkit didn't make them in the first place, of course...
              string smi = DACLIB::create_cansmi( *prod_mol );
              if( verbose ) {
                cout << "AWOOGA - got some radicals for " << in_mol.GetTitle() << " : " << smi << endl;
              }
              delete prod_mol;
            } else {
              // fix any chiral centres that may have been affected by reaction
              remove_altered_stereochem( libgen , prod_mol );
              string smi = DACLIB::create_cansmi( *prod_mol );
              if( all_can_smis.find( smi ) == all_can_smis.end() ) {
                if( add_smirks_to_name ) {
                  string curr_name = prod_mol->GetTitle();
                  curr_name += string( " " ) + smirks_[smirks_num].first;
                  prod_mol->SetTitle( curr_name );
                }
                ret_mols.push_back( prod_mol );
                if( ret_mols.size() > max_out_mols_ ) {
                  // it's going to take too long
                  for( size_t j = 0 , js = ret_mols.size() ; j < js ; ++j ) {
                    delete ret_mols[j];
                  }
                  ret_mols.clear();
                  throw TooManyOutMols( in_mol );
                }
                all_can_smis.insert( smi );
                if( verbose ) {
                  cout << endl << "New product in tautomer enumerator : " << smi << endl
                       << "Made from " << DACLIB::create_cansmi( *ret_mols[i] ) << endl
                       << "Using SMIRKS : " << smirks_[smirks_num].first << " : " << smirks_[smirks_num].second << endl
                       << "Expanded to : " << exp_smirks_[smirks_num] << endl;
                }
              } else {
                delete prod_mol; // we've already got this molecule
              }
            }
          }
        } else {
#ifdef NOTYET
          cout << "No prods for this libgen" << endl;
#endif
        }
        ++smirks_num; // counter for libgen/smirks, for debugging
      }
    }
#ifdef NOTYET
    cout << "Number of tautomers currently : " << ret_mols.size() << endl;
    BOOST_FOREACH( OEMolBase *rm , ret_mols ) {
      cout << DACLIB::create_cansmi( *rm ) << endl;
    }
#endif

    if( start_size == ret_mols.size() ) {
      break;
    } else {
      next_start = start_size;
    }
  }

  // put molecules in consistent order
  vector<pair<string,OEMolBase *> > smiles;
  create_smiles( ret_mols , smiles );

  ret_mols.clear();
  transform( smiles.begin() , smiles.end() ,
             back_inserter( ret_mols ) ,
             bind( &pair<string,OEMolBase *>::second, _1 ) );

  return ret_mols;

}

// ****************************************************************************
vector<string> TautEnum::enumerate_smiles( OEMolBase &in_mol , bool verbose ,
                                           bool add_smirks_to_name ) {

  vector<OEMolBase *> all_mols = enumerate( in_mol , verbose , add_smirks_to_name );

  vector<pair<string,OEMolBase *> > smiles;
  create_smiles( all_mols , smiles );

  vector<string> ret_val;
  transform( smiles.begin() , smiles.end() ,
             back_inserter( ret_val ) ,
             bind( &pair<string,OEMolBase *>::first , _1 ) );

  for( size_t i = 1 , is = smiles.size() ; i < is ; ++i ) {
    delete smiles[i].second;
  }

  return ret_val;

}

// ************************************************************************************
// remove any stereochemistry from atoms affected by the reaction
void TautEnum::remove_altered_stereochem( pOELibGen &libgen , OEMolBase *mol ) {

  for( int i = 0 , is = libgen->NumReactants() ; i < is ; ++i ) {
    OEIter<OEMolBase> sms = libgen->GetStartingMaterial( i );
    for( ; sms ; ++sms ) {
      // fix atom stereo
      for( OEIter<OEAtomBase> atom = sms->GetAtoms( OEHasMapIdx() ) ; atom ; ++atom ) {
        if( atom->HasStereoSpecified() ) {
          // cout << "stereo on atom " << atom->GetIdx() << " : " << atom->GetMapIdx() << " for " << sms->GetTitle() << endl;
          OEIter<OEAtomBase> patom = mol->GetAtoms( OEHasMapIdx( atom->GetMapIdx() ) );
          if( patom->HasStereoSpecified( OEAtomStereo::Tetra ) ) {
            if( atom->GetAtomicNum() != patom->GetAtomicNum() ||
                atom->GetDegree() != patom->GetDegree() ||
                atom->GetHvyDegree() != patom->GetHvyDegree() ||
                atom->GetValence() != patom->GetValence() ||
                atom->GetHvyDegree() != patom->GetHvyDegree() ||
                atom->GetHyb() != patom->GetHyb() ||
                atom->GetTotalHCount() != patom->GetTotalHCount() ) {
              patom->SetStereo( vector<OEAtomBase *>() , OEAtomStereo::Tetra , OEAtomStereo::Undefined );
              cout << "Stereochem removed for tautomer of " << sms->GetTitle() << endl;
            }
          }
        }
      }
#ifdef NOTYET
      // fix bond stereo - doesn't seem to be required, but might as well leave the
      // testing code that established this, just in case.
      for( OEIter<OEBondBase> bond = sms->GetBonds() ; bond ; ++bond ) {
        if( bond->GetBgn()->GetMapIdx() && bond->GetEnd()->GetMapIdx() ) {
          cout << "Bond between mapped atoms " << bond->GetBgn()->GetIdx() << " , " << bond->GetBgn()->GetMapIdx()
               << " and " << bond->GetEnd()->GetIdx() << " , " << bond->GetEnd()->GetMapIdx() << endl;
          cout << "Stereo specified : " << bond->HasStereoSpecified( OEBondStereo::CisTrans ) << endl;
          OEAtomBase *pab = mol->GetAtom( OEHasMapIdx( bond->GetBgn()->GetMapIdx() ) );
          OEAtomBase *pae = mol->GetAtom( OEHasMapIdx( bond->GetEnd()->GetMapIdx() ) );
          if( pab && pae ) {
            OEBondBase *pbond = mol->GetBond( pab , pae );
            if( pbond ) {
              cout << "P bond stereo : " << pbond->HasStereoSpecified( OEBondStereo::CisTrans ) << endl;
            } else {
              cout << "No product bond" << endl;
            }
          }
        }
      }
#endif
    }
  }

}

// ****************************************************************************
void create_smiles( vector<OEMolBase *> &all_mols ,
                    vector<pair<string,OEMolBase *> > &smiles ) {

  BOOST_FOREACH( OEMolBase *mol , all_mols ) {
    string smi;
    OECreateSmiString( smi , *mol , OESMILESFlag::Canonical | OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo );
    smiles.push_back( make_pair( smi , mol ) );
  }

  sort( smiles.begin() , smiles.end() ,
        bind( greater<string>() ,
              bind( &pair<string,OEMolBase *>::first , _1 ) ,
              bind( &pair<string,OEMolBase *>::first , _2 ) ) );

}
