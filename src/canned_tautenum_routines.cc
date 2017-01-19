//
// file canned_tautenum_routines.cc
// David Cosgrove
// AstraZeneca
// 16th February 2012
//
// This file contains a few functions for doing standard tautomer transformations
// in one call. They are intended to be used in other programs either linked directly
// in C++ or as a Python module.

#include "TautEnum.H"
#include "TautStand.H"
#include "taut_enum_default_vector_bindings.H"
#include "taut_enum_default_standardise_smirks.H"
#include "taut_enum_default_enum_smirks_extended.H"
#include "taut_enum_default_enum_smirks_orig.H"
#include "taut_enum_protonate_a.H"
#include "taut_enum_protonate_b.H"
#include "taut_enum_protonate_vb.H"

#include <oechem.h>

#include <iostream>
#include <vector>

using namespace std;
using namespace OEChem;

// ****************************************************************************
void prepare_molecule( OEMolBase &mol ) {

  OEPerceiveChiral( mol );
  OEAssignAromaticFlags( mol );
  // put the input molecule into a canonical form, so that transformations
  // for which the output is ambiguous (e.g. CCC(=N)C -> CC=C(N)C or CCC(N)=C
  // come out consistently
  string can_smi , mol_name;
  mol_name = mol.GetTitle();
  OECreateSmiString( can_smi , mol , OESMILESFlag::Canonical | OESMILESFlag::AtomStereo | OESMILESFlag::BondStereo );
  mol.Clear();
  OEParseSmiles( mol , can_smi );
  mol.SetTitle( mol_name );

}


// *********************************************************************************
// Take the molecule and produce a standardised tautomer, suitable for input into
// the other routines
OEMolBase *standardise_tautomer( OEMolBase &in_mol ) {

  static TautStand *taut_stand = 0;
  if( !taut_stand ) {
    taut_stand = new TautStand( DACLIB::STAND_SMIRKS , DACLIB::VBS );
  }

  prepare_molecule( in_mol );
  OEMolBase *ret_mol = taut_stand->standardise( in_mol , false );

  return ret_mol;

}

// *********************************************************************************
vector<OEMolBase *> enumerate_ions( OEMolBase &in_mol ,
                                    const string &prot_stand_smirks ,
                                    const string &prot_enum_smirks ,
                                    const string &prot_smirks_vbs ) {

  static TautStand *prot_stand = 0;
  if( !prot_stand ) {
    if( prot_stand_smirks.empty() ) {
      prot_stand = new TautStand( DACLIB::PROTONATE_A , DACLIB::SET_PROT_VB );
    } else {
      prot_stand = new TautStand( prot_stand_smirks , prot_smirks_vbs );
    }
  }
  static TautEnum *prot_enum = 0;
  if( !prot_enum ) {
    if( prot_stand_smirks.empty() ) {
      prot_enum = new TautEnum( DACLIB::PROTONATE_B , DACLIB::SET_PROT_VB );
    } else {
      prot_enum = new TautEnum( prot_enum_smirks , prot_smirks_vbs );
    }
  }

  OEMolBase *std_mol = prot_stand->standardise( in_mol , false , false , false );

  vector<OEMolBase *> ret_mols;
  try {
    ret_mols = prot_enum->enumerate( *std_mol , false , false );
  } catch( TooManyOutMols &e ) {
    cerr << "Maximum number of tautomers generated for " << in_mol.GetTitle()
         << " so none generated." << endl;
    ret_mols.push_back( OENewMolBase( *std_mol , OEMolBaseType::OEDefault ) );
    return ret_mols;
  }

  delete std_mol;
  return ret_mols;

}

// *********************************************************************************
vector<OEMolBase *> enumerate_tautomers( OEMolBase &in_mol ,
                                         const string &smirks_defs ,
                                         const string &smirks_vbs ) {

  static TautEnum *taut_enum = 0;
  if( !taut_enum ) {
    if( smirks_defs.empty() ) {
      taut_enum = new TautEnum( DACLIB::ENUM_SMIRKS_EXTENDED , DACLIB::VBS );
    } else {
      taut_enum = new TautEnum( smirks_defs , smirks_vbs );
    }
  }

  OEMolBase *std_mol = standardise_tautomer( in_mol );
  vector<OEMolBase *> taut_mols;
  try {
    // false for not verbose output
    taut_mols = taut_enum->enumerate( *std_mol , false );
  } catch( TooManyOutMols &e ) {
    cerr << "Maximum number of tautomers generated for " << in_mol.GetTitle() << " so none generated." << endl;
    taut_mols = vector<OEMolBase *>( 1 , OENewMolBase( in_mol , OEMolBaseType::OEDefault ) );
  }
  delete std_mol;

  return taut_mols;

}

// *********************************************************************************
OEMolBase *canonical_tautomer( OEMolBase &in_mol ) {

  vector<OEMolBase *> all_tauts = enumerate_tautomers( in_mol , string( "" ) ,
						       string( "" ) );
  vector<pair<string,OEMolBase *> > smiles;

  // smiles and all_tauts will come out of this pointing at the same OEMolBases
  create_smiles( all_tauts , smiles );

  // delete the pointer from smiles rather than all_tauts because they are in different
  // orders smiles is the one that counts.
  for( size_t i = 1 , is = smiles.size() ; i < is ; ++i ) {
    delete smiles[i].second;
  }

  return smiles[0].second;

}

// *********************************************************************************
// generate set of SMILES strings for tautomers of in_smi. First entry in
// return vector will be a canonical SMILES for the input string, followed
// by canonical SMILES of all other tautomers.
vector<string> enumerate_tautomers_smiles( const string &in_smi ) {

  OEGraphMol in_mol;
  OEParseSmiles( in_mol , in_smi );
  vector<string> ret_val;

  unsigned int oeflavour = OESMILESFlag::ISOMERIC ^ OESMILESFlag::AtomMaps;
  string smi;
  OECreateSmiString( smi , in_mol , oeflavour );
  ret_val.push_back( smi );

  vector<OEMolBase *> tauts = enumerate_tautomers( in_mol , string( "" ) ,
						   string( "" ) );
  for( size_t i = 0 , is = tauts.size() ; i < is ; ++i ) {
    OECreateSmiString( smi , *tauts[i] , oeflavour );
    if( ret_val.end() == find( ret_val.begin() , ret_val.end() , smi ) ) {
      ret_val.push_back( smi );
    }
    delete tauts[i];
  }

  return ret_val;

}
