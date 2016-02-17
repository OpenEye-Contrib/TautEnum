//
// File TautEnumCallable.cc
// David Cosgrove
// AstraZeneca
// 8th February 2012.
//

#include "TautEnum.H"
#include "TautStand.H"
#include "TautEnumCallableBase.H"
#include "FileExceptions.H"
#include "taut_enum_default_vector_bindings.H"
#include "taut_enum_default_standardise_smirks.H"
#include "taut_enum_default_enum_smirks_extended.H"
#include "taut_enum_default_enum_smirks_orig.H"
#include "taut_enum_protonate_a.H"
#include "taut_enum_protonate_b.H"
#include "taut_enum_protonate_vb.H"

#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/lexical_cast.hpp>

#include <oechem.h>

#include <vector>

using namespace OEChem;
using namespace std;

namespace DACLIB {
void apply_daylight_aromatic_model( OEMolBase &mol );
string create_cansmi( const OEMolBase &in_mol );
}

// in canned_tautenum_routines.cc
void prepare_molecule( OEMolBase &mol );

// ****************************************************************************
void TautEnumCallableBase::operator ()() {

  TautStand *taut_stand = 0;
  TautEnum *taut_enum = 0;

  // create this standardise/enumerate pair as we're always going to standardise.
  const string &enum_smirks = tes_.extended_enumeration() ? DACLIB::ENUM_SMIRKS_EXTENDED : DACLIB::ENUM_SMIRKS_ORIG;
  create_enumerator_objects( tes_.standardise_smirks_file() , tes_.enumerate_smirks_file() ,
                             tes_.vb_file() , DACLIB::STAND_SMIRKS , enum_smirks ,
                             DACLIB::VBS , taut_stand , taut_enum );

  TautStand *prot_stand = 0;
  TautEnum *prot_enum = 0;
  if( tes_.enumerate_protonation() ) {
    create_enumerator_objects( tes_.protonation_standardisation_file() , tes_.protonation_enumeration_file() ,
                               tes_.protonation_vb_file() , DACLIB::PROTONATE_A , DACLIB::PROTONATE_B ,
                               DACLIB::SET_PROT_VB , prot_stand , prot_enum );
  }

  OEMolBase *in_mol = OENewMolBase( OEMolBaseType::OEDefault );
  int mol_num = 0;

  while( read_next_molecule( *in_mol ) ) {
    ++mol_num;
    if( tes_.verbose() ) {
      cout << "Processing " << in_mol->GetTitle() << " : " << DACLIB::create_cansmi( *in_mol ) << " (" << mol_num << ")"  << endl;
    }
    // OEReadMolecule doesn't do quite as much of a setup of the molecules,
    // as I recall. Do it explicitly, just to be safer.
    prepare_molecule( *in_mol );
    if( tes_.verbose() ) {
      cout << "Pre-processed molecule : " <<  DACLIB::create_cansmi( *in_mol ) << endl;
    }

    vector<OEMolBase *> out_mols;
    OEMolBase *std_mol = 0;
    if( taut_stand ) {
      std_mol = taut_stand->standardise( *in_mol , tes_.verbose() ,
                                         tes_.add_smirks_to_name() ,
                                         tes_.strip_salts() );
    } else {
      std_mol = OENewMolBase( *in_mol , OEMolBaseType::OEDefault );
    }

    if( !tes_.standardise_only() ) {
      if( tes_.extended_enumeration() || tes_.original_enumeration() ) {
        try {
          vector<OEMolBase *> taut_mols = taut_enum->enumerate( *std_mol , tes_.verbose() ,
                                                                tes_.add_smirks_to_name() );
          out_mols.insert( out_mols.end() , taut_mols.begin() , taut_mols.end() );
        } catch( TooManyOutMols &e ) {
          // just leave it as the standardised molecule
          cerr << "Maximum number of tautomers generated for " << in_mol->GetTitle() << " so none generated." << endl;
          out_mols.push_back( OENewMolBase( *std_mol , OEMolBaseType::OEDefault ) );
          if( tes_.add_smirks_to_name() ) {
            string new_name = in_mol->GetTitle() + string( " __MAX_TAUTS__" );
            out_mols.back()->SetTitle( new_name );
          }
        }
      }

      if( prot_enum ) {
        if( out_mols.empty() ) {
          // just doing an enumerate_protonation job. May need to do strip salts.
          OEMolBase *std_prot_mol = prot_stand->standardise( *std_mol , tes_.verbose() ,
                                                             tes_.add_smirks_to_name() ,
                                                             true );
          try {
            vector<OEMolBase *> prot_mols = prot_enum->enumerate( *std_prot_mol , tes_.verbose() ,
                                                                  tes_.add_smirks_to_name() );
            out_mols.insert( out_mols.end() , prot_mols.begin() , prot_mols.end() );
          } catch( TooManyOutMols &e ) {
            cerr << "Maximum number of ionisation states generated for " << in_mol->GetTitle() << " so none generated." << endl;
            // just leave it as it was. I think it's pretty unlikely to happen.
            std_prot_mol = 0;
          }
          delete std_prot_mol;
        } else {
          // in this case, we don't want to include the output from the tautomer enumeration
          // in the output, but we do want to pass each tautomer through the protonation
          // enumerator
          vector<OEMolBase *> prot_out_mols;
          for( int i = 0 , is = out_mols.size() ; i < is ; ++i ) {
            // strip_salts will already have been applied by taut_stand if we wanted to do it,
            // as all input mols are standardised.
            OEMolBase *std_prot_mol = prot_stand->standardise( *out_mols[i] , tes_.verbose() ,
                                                               tes_.add_smirks_to_name() ,
                                                               false );
            try {
              vector<OEMolBase *> prot_mols = prot_enum->enumerate( *std_prot_mol , tes_.verbose() ,
                                                                    tes_.add_smirks_to_name() );
              prot_out_mols.insert( prot_out_mols.end() , prot_mols.begin() , prot_mols.end() );
            } catch( TooManyOutMols &e ) {
              cerr << "Maximum number of ionisation states generated for " << in_mol->GetTitle() << " tautomer " << i << " so none generated." << endl;
              // just leave it as it was. I think it's pretty unlikely to happen.
            }
            delete std_prot_mol;
          }
          // empty out_mols and replace with prot_out_mols
          for( int i = 0 , is = out_mols.size() ; i < is ; ++i ) {
            delete out_mols[i];
          }
          out_mols = prot_out_mols;
        }
      }
    } else {
#ifdef NOTYET
      cout << "standardise only" << endl;
#endif
      out_mols.push_back( OENewMolBase( *std_mol , OEMolBaseType::OEDefault ) );
    }

    sort_and_uniquify_molecules( out_mols );

    if( tes_.include_input_in_output() ) {
      write_molecule( *in_mol );
    }
    if( tes_.canonical_tautomer() ) {
      if( out_mols.empty() ) {
        // probably hit the exception for too many tautomers, so write standardised input mol
        out_mols.push_back( OENewMolBase( *std_mol , OEMolBaseType::OEDefault ) );
      }
      if( tes_.add_numbers_to_name() ) {
        // it's not very sensible, but the user might ask for it
        out_mols.front()->SetTitle( out_mols.front()->GetTitle() + tes_.name_postfix() + string( "1" ) );
      }
      write_molecule( *out_mols.front() );
    } else {
      output_molecules( out_mols );
    }
    while( !out_mols.empty() ) {
#ifdef NOTYET
      cout << "deleting out_mols : " << out_mols.size() << endl;
#endif
      delete out_mols.back();
      out_mols.pop_back();
    }
    in_mol->Clear();
    delete std_mol;
  }

  // give it time to clear up properly
  boost::this_thread::sleep( boost::posix_time::seconds( 1 ) );
  delete in_mol;
  delete taut_stand;
  delete taut_enum;
  delete prot_stand;
  delete prot_enum;

}

// ****************************************************************************
// make the TautStand and TautEnum objects, using the relevant data from tes_
void TautEnumCallableBase::create_enumerator_objects( const string &stand_smirks_file ,
                                                      const string &enum_smirks_file ,
                                                      const string &vb_file ,
                                                      const string &default_stand_smirks ,
                                                      const string &default_enum_smirks ,
                                                      const string &default_vbs ,
                                                      TautStand *&taut_stand ,
                                                      TautEnum *&taut_enum ) {


  if( !stand_smirks_file.empty() ) {
#ifdef NOTYET
    cout << "standardise SMIRKS file : " << stand_smirks_file << endl;
#endif
    try {
      // true just forces the file reading c'tor to be used.
      taut_stand = new TautStand( stand_smirks_file , vb_file , true );
    } catch( DACLIB::FileReadOpenError &e ) {
      cout << e.what() << endl;
      exit( 1 );
    }
  } else {
#ifdef NOTYET
    cout << "default standardisation SMIRKS" << endl;
#endif
    taut_stand = new TautStand( default_stand_smirks , default_vbs ); // using default standardisation SMIRKS
  }

  if( !enum_smirks_file.empty() ) {
#ifdef NOTYET
    cout << "enumerate SMIRKS file : " << enum_smirks_file << endl;
#endif
    try {
      taut_enum = new TautEnum( enum_smirks_file , vb_file , true ,
                                tes_.max_tautomers() );
    } catch( DACLIB::FileReadOpenError &e ) {
      cout << e.what() << endl;
      exit( 1 );
    }
  } else {
#ifdef NOTYET
    cout << "default enumeration SMIRKS" << endl;
#endif
    taut_enum = new TautEnum( default_enum_smirks , default_vbs , tes_.max_tautomers() );
  }

}

// ****************************************************************************
void TautEnumCallableBase::output_molecules( vector<OEMolBase *> &out_mols ) {

  for( int i = 0 , is = out_mols.size() ; i < is ; ++i ) {
    if( tes_.add_numbers_to_name() ) {
      out_mols[i]->SetTitle( out_mols[i]->GetTitle() + tes_.name_postfix() + boost::lexical_cast<string>( i + 1 ) );
    }
    write_molecule( *out_mols[i] );
  }

}

// ****************************************************************************
void sort_and_uniquify_molecules( vector<OEMolBase *> &mols ) {

  vector<pair<string,OEMolBase *> > tmp_mols;
  create_smiles( mols , tmp_mols );

  for( int i = 1 , is = tmp_mols.size() ; i < is ; ++i ) {
    if( tmp_mols[i].first == tmp_mols[i-1].first ) {
      delete tmp_mols[i].second;
      tmp_mols[i].second = static_cast<OEMolBase *>( 0 );
    }
  }

  using namespace boost;
  tmp_mols.erase( remove_if( tmp_mols.begin() , tmp_mols.end() ,
                             bind( logical_not<bool>() ,
                                   bind( &pair<string,OEMolBase *>::second , _1 ) ) ) ,
                  tmp_mols.end() );

  mols.clear();
  transform( tmp_mols.begin() , tmp_mols.end() , back_inserter( mols ) ,
             bind( &pair<string,OEMolBase *>::second , _1 ) );

}
