//
// file DACOEMolAtomIndex.H
// David Cosgrove
// AstraZeneca
// 3rd March 2008
//
// The number returned by OEAtomBase->GetIdx() (OE Atom Idx) is not preserved
// when the molecule is copied. I got so in the habit of believing that it was,
// and relied on this for so many things, that I have implemented something
// that is. If the indices are accessed without being created, they will be
// created on the fly.
// Max atom index is intended always to be 1 more than the highest index
// in the molecule, so it can be used sensibly for sizing arrays.

#include <oechem.h>

// **********************************************************************
namespace DACLIB {

  static const int ATOM_INDEX_TAG = OESystem::OEGetTag( "ATOM_INDEX_TAG" );
  static const int MAX_ATOM_INDEX_TAG = OESystem::OEGetTag( "MAX_ATOM_INDEX_TAG" );

  static unsigned int max_atom_index( OEChem::OEMolBase &mol );
  static void set_atom_index( OEChem::OEAtomBase &atom ,
				  unsigned int new_ind );
  static void create_atom_indices( OEChem::OEMolBase &mol );
  static unsigned int atom_index( const OEChem::OEAtomBase &atom );

  // *******************************************************************
  // Creates the atom index if it hasn't been done yet
  static unsigned int max_atom_index( OEChem::OEMolBase &mol ) {

    if( !mol.HasData( MAX_ATOM_INDEX_TAG ) ) {
      create_atom_indices( mol );
    }

    return mol.GetData<unsigned int>( MAX_ATOM_INDEX_TAG );

  }

  // *******************************************************************
  // Set the atom index to the given value, updating max atom index for
  // the parent molecule if necessary.
  static void set_atom_index( OEChem::OEAtomBase &atom ,
                              unsigned int new_ind ) {

    atom.DeleteData( ATOM_INDEX_TAG ); // in case we're re-indexing
    atom.SetData<unsigned int>( ATOM_INDEX_TAG , new_ind );
    if( new_ind >= max_atom_index( *atom.GetParent() ) ) {
      atom.GetParent()->SetData( MAX_ATOM_INDEX_TAG , new_ind + 1 );
    }

  }

  // *******************************************************************
  static void create_atom_indices( OEChem::OEMolBase &mol ) {

    mol.SetData<unsigned int>( MAX_ATOM_INDEX_TAG , mol.GetMaxAtomIdx() + 1 );
    for( OESystem::OEIter<OEChem::OEAtomBase> atom =
	   mol.GetAtoms() ; atom ; ++atom ) {
      set_atom_index( atom , atom->GetIdx() );
    }

  }

  // *******************************************************************
  // Creates the atom index if it hasn't been done yet. The const-ness of
  // the reference doesn't seem to matter even though the state of the
  // molecule can be changed.
  static unsigned int atom_index( const OEChem::OEAtomBase &atom ) {

    if( !atom.HasData( ATOM_INDEX_TAG ) ) {
      create_atom_indices( *atom.GetParent() );
    }
    
    return atom.GetData<unsigned int>( ATOM_INDEX_TAG );
    
  }

  // *******************************************************************
  // Unary predicate to identify atom with given index tag - borrowed
  // heavily from code from Dave Riley
  class HasAtomIndex : public OESystem::OEUnaryPredicate<OEChem::OEAtomBase> {

  public :

    HasAtomIndex( unsigned int id ) : search_id_( id ) {}

    // we'll be lazy, and rely on the default copy c'tor and operator=.
    // DR wouldn't be, but for now a bitwise copy is fine.

    // need to make a virtual function CreateCopy.
    virtual OESystem::OEUnaryPredicate<OEChem::OEAtomBase> *CreateCopy() const {
      return new HasAtomIndex( *this );
    }

    // and what we're really after
    virtual bool operator()( const OEChem::OEAtomBase &atom ) const {
      return atom_index( atom ) == search_id_;
    }

  private :

    unsigned int search_id_; // set in c'tor

  };

} // EO namespace DACLIB
