//
// file apply_daylight_arom_model_to_oemol.cc
// David Cosgrove
// AstraZeneca
// 27th January 2009
//

#include <oechem.h>

using namespace OEChem;
using namespace OESystem;

namespace DACLIB {

// ****************************************************************************
void apply_daylight_aromatic_model( OEMolBase &mol ) {

  OESuppressHydrogens( mol );
  OEClearAromaticFlags( mol );

  // we need to set the integer bond types correctly before calling OEKekulize,
  // and we call OEKekulise to make sure that the aromatic flags are
  // truly hosed.  I assume that at some point I gained empirical evidence
  // that OEClearAromaticFlags wasn't strong enough voodoo.
  // mols read from most file types will not need this, but OEB for one does
  // so do it for all. It's not going to be the rate-limiting step.
  for( OEIter<OEBondBase> bond = mol.GetBonds(); bond ; ++bond ) {
    if( bond->IsAromatic() ) {
      bond->SetIntType( 5 );
    } else {
      bond->SetIntType( bond->GetOrder() );
    }
  }
  OEKekulize( mol );
  OEAssignAromaticFlags( mol , OEAroModelDaylight , true );
  OESuppressHydrogens( mol ); // we really want those hydrogens suppressed!

}

} // EO namespace DACLIB
