//
// file DACOEMolBondIndex.H
// David Cosgrove
// AstraZeneca
// 3rd March 2008
//
// The number returned by OEBondBase->GetIdx() (OE Bond Idx) is not preserved
// when the molecule is copied. I got so in the habit of believing that it was,
// and relied on this for so many things, that I have implemented something
// that is. If the indices are accessed without being created, they will be
// created on the fly.
// Max Bond index is intended always to be 1 more than the highest index
// in the molecule, so it can be used sensibly for sizing arrays.

#include <oechem.h>

// **********************************************************************
namespace DACLIB {

  static const int BOND_INDEX_TAG = OESystem::OEGetTag( "BOND_INDEX_TAG" );
  static const int MAX_BOND_INDEX_TAG = OESystem::OEGetTag( "MAX_BOND_INDEX_TAG" );

  static unsigned int max_bond_index( OEChem::OEMolBase &mol );
  static void set_bond_index( OEChem::OEBondBase &bond ,
                              unsigned int new_ind );
  static void create_bond_indices( OEChem::OEMolBase &mol );
  static unsigned int bond_index( const OEChem::OEBondBase &bond );

  // *******************************************************************
  // Creates the bond index if it hasn't been done yet
  static unsigned int max_bond_index( OEChem::OEMolBase &mol ) {

    if( !mol.HasData( MAX_BOND_INDEX_TAG ) ) {
      create_bond_indices( mol );
    }

    return mol.GetData<unsigned int>( MAX_BOND_INDEX_TAG );

  }

  // *******************************************************************
  // Set the bond index to the given value, updating max bond index for
  // the parent molecule if necessary.
  static void set_bond_index( OEChem::OEBondBase &bond ,
                              unsigned int new_ind ) {

    bond.SetData<unsigned int>( BOND_INDEX_TAG , new_ind );
    if( new_ind >= max_bond_index( *bond.GetParent() ) ) {
      bond.GetParent()->SetData( MAX_BOND_INDEX_TAG , new_ind + 1 );
    }

  }

  // *******************************************************************
  static void create_bond_indices( OEChem::OEMolBase &mol ) {

    mol.SetData<unsigned int>( MAX_BOND_INDEX_TAG ,
                   mol.GetMaxBondIdx() + 1 );
    for( OESystem::OEIter<OEChem::OEBondBase> bond =
       mol.GetBonds() ; bond ; ++bond ) {
      set_bond_index( bond , bond->GetIdx() );
    }

  }

  // *******************************************************************
  // Creates the bond index if it hasn't been done yet. The const-ness of
  // the reference doesn't seem to matter even though the state of the
  // molecule can be changed.
  static unsigned int bond_index( const OEChem::OEBondBase &bond ) {

    if( !bond.HasData( BOND_INDEX_TAG ) ) {
      create_bond_indices( *bond.GetParent() );
    }
    
    return bond.GetData<unsigned int>( BOND_INDEX_TAG );
    
  }

  // *******************************************************************
  // Unary predicate to identify bond with given index tag - borrowed
  // heavily from code from Dave Riley
  class HasBondIndex : public OESystem::OEUnaryPredicate<OEChem::OEBondBase> {

  public :

    HasBondIndex( unsigned int id ) : search_id_( id ) {}

    // we'll be lazy, and rely on the default copy c'tor and operator=.
    // DR wouldn't be, but for now a bitwise copy is fine.

    // need to make a virtual function CreateCopy.
    virtual OESystem::OEUnaryPredicate<OEChem::OEBondBase> *CreateCopy() const {
      return new HasBondIndex( *this );
    }

    // and what we're really after
    virtual bool operator()( const OEChem::OEBondBase &bond ) const {
      return bond_index( bond ) == search_id_;
    }

  private :

    unsigned int search_id_; // set in c'tor

  };

} // EO namespace DACLIB
