//
// file TautStand.H
// David Cosgrove
// AstraZeneca
// 2nd August 2011
//
// This class takes two file names, of SMIRKS and vector bindings to be
// used for standardising the tautomer of a molecule. The SMIRKS are
// expected to be derived from Pete Kenny's original Leatherface
// work, but any old SMIRKS will work.  They are applied in the
// order they are in the file, each one being applied exhaustively
// to the molecule before moving on. There's therefore a clear
// dependence on the order of the SMIRKS in the input file.

#ifndef TAUTSTAND_H
#define TAUTSTAND_H

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

// ****************************************************************************

namespace OEChem {
class OEMolBase;
class OELibraryGen;
}

typedef boost::shared_ptr<OEChem::OELibraryGen> pOELibGen;

// ****************************************************************************

class TautStand {

public :

  // load from strings, probably from a default include file
  TautStand( const std::string &smirks_string , const std::string &vb_string );
  // load from SMILES, with a dummy bool to make the signatures different
  TautStand( const std::string &smirks_file , const std::string &vb_file,
             bool dummy );
  TautStand( const TautStand &rhs ); // needs copy c'tor for threading

  OEChem::OEMolBase *standardise( OEChem::OEMolBase &in_mol , bool verbose = false ,
                                  bool add_smirks_to_name = false ,
                                  bool strip_salts = false );

private :

  std::string smirks_file_;
  std::string vb_file_;
  std::vector<std::pair<std::string,std::string> > smirks_;
  std::vector<std::pair<std::string,std::string> > vbs_;
  std::vector<std::string> exp_smirks_; // the SMIRKS with vector bindings expanded
  std::vector<pOELibGen> lib_gens_; // the reaction objects, built from the SMIRKS

};

#endif // TAUTSTAND_H
