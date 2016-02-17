//
// file radical_atoms.cc
// David Cosgrove
// AstraZeneca
// 15th February 2013
//
// This function takes a molecule and builds a list of those atoms that
// are free radicals, taken from a selected list of elements.

#include <map>
#include <vector>

#include <oechem.h>
#include <oedepict.h>

#include <boost/tuple/tuple.hpp>

using namespace OEChem;
using namespace OEDepict;
using namespace OESystem;

namespace DACLIB {

// *****************************************************************************
void radical_atoms( OEMolBase &mol ,
                    std::vector<OEAtomBase *> &rad_atoms ) {

  static std::map<unsigned int,boost::tuple<int,int,int> > ELECTRONS;
  if( ELECTRONS.empty() ) {
    ELECTRONS[OEElemNo::C] = boost::tuple<int,int,int>( 4 , -1 , -1 );
    ELECTRONS[OEElemNo::N] = boost::tuple<int,int,int>( 5 , -1 , -1 );
    ELECTRONS[OEElemNo::O] = boost::tuple<int,int,int>( 6 , -1 , -1 );
    ELECTRONS[OEElemNo::Si] = boost::tuple<int,int,int>( 4 , -1 , -1 );
    ELECTRONS[OEElemNo::P] = boost::tuple<int,int,int>( 5 , 3 , -1 );
    ELECTRONS[OEElemNo::S] = boost::tuple<int,int,int>( 6 , 4 , 2 );
  }

  for( OEIter<OEAtomBase> atom = mol.GetAtoms() ; atom ; ++atom ) {
    std::map<unsigned int,boost::tuple<int,int,int> >::iterator p = ELECTRONS.find( atom->GetAtomicNum() );
    if( p != ELECTRONS.end() ) {
      bool rad = true;
      int n = p->second.get<0>() + atom->GetValence() - atom->GetFormalCharge();
      if( !(n - 8) ) {
        rad = false;
      }
      if( -1 != p->second.get<1>() ) {
        n = p->second.get<1>() + atom->GetValence() - atom->GetFormalCharge();
        if( !(n - 8) ) {
          rad = false;
        }
      }
      if( -1 != p->second.get<2>() ) {
        n = p->second.get<2>() + atom->GetValence() - atom->GetFormalCharge();
        if( !(n - 8) ) {
          rad = false;
        }
      }
      if( rad ) {
        rad_atoms.push_back( atom );
      }
    }
  }

}

#ifdef NOTYET
// *****************************************************************************
// this one is just for checking the previous one, really.  It is almost
// certainly slower as it has to do a full depiction calculation.
// It's not compiled in by default as it requires an extra library (oedepict)
// that isn't always linked to non-graphical programs.
bool radical_atoms( OEMolBase &mol ) {

  OEPrepareDepiction(mol);

  unsigned int width  = 200;
  unsigned int height = 200;
  OE2DMolDisplayOptions opts(width, height, OEScale::AutoScale);

  OE2DMolDisplay disp(mol, opts);

  for (OEIter<OE2DAtomDisplay> adisp = disp.GetAtomDisplays(); adisp; ++adisp) {
    if( adisp->GetRadical() ) {
      return true;
    }
  }

  return false;

}
#endif

} // EO namespace DACLIB
