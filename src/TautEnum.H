//
// file TautEnum.H
// David Cosgrove
// AstraZeneca
// 3rd August 2011
//
// This class does a tautomer enumeration on a molecule.
// It takes two file names, a SMIRKS file and an associated file
// of vector bindings.  The SMIRKS are expected to be derived from
// Pete Kenny's original Leatherface work, but any old SMIRKS should
// work.  For optimum performance, the molecule should be in a
// standard tautomer, as produced by the companion class TautStand,
// for example.
// Following Pete's work, a canonical tautomer can be returned. This
// is no more than the first tautomer when the canonical SMILES for
// all tautomers are sorted in ascending alphanumeric order.
//
// It works best if the input molecule has been standardised using
// the associacted class TautStand.

#ifndef TAUTENUM_H
#define TAUTENUM_H

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

class TautEnum {

public :

  TautEnum( const std::string &smirks_string , const std::string &vbs_string ,
            unsigned int max_t = 256 ); // load the canned defaults
  // dummy to distinguish signatures
  TautEnum( const std::string &smirks_file , const std::string &vb_file ,
            bool dummy , unsigned int max_t = 256 );
  TautEnum( const TautEnum &rhs ); // needs copy c'tor for threading

  // tautomers() throws an exception of type TooManyTautomers if max_tauts_ is exceeded
  // and the other 3 functions below call tautomers don't catch the exception.
  // tautomers always has the input tautomer first in the output vector.
  // likewise tautomer_smiles.
  std::vector<OEChem::OEMolBase *> enumerate( OEChem::OEMolBase &in_mol , bool verbose = false ,
                                              bool add_smirks_to_name = false );
  std::vector<std::string> enumerate_smiles( OEChem::OEMolBase &in_mol , bool verbose = false ,
                                             bool add_smirks_to_name = false );

  unsigned int max_out_mols() const { return max_out_mols_; }

private :

  std::string smirks_file_;
  std::string vb_file_;
  std::vector<std::pair<std::string,std::string> > smirks_;
  std::vector<std::pair<std::string,std::string> > vbs_;
  std::vector<std::string> exp_smirks_; // with vector bindings expanded
  std::vector<pOELibGen> lib_gens_; // the SMIRKS transformed into reaction objects
  const unsigned int max_out_mols_; //  maximum number of tautomers to be generated. Returns just the input molecule (i.e. no tautomers) if exceeded.

  // remove any stereochemistry from atoms affected by the reaction
  void remove_altered_stereochem( pOELibGen &libgen , OEChem::OEMolBase *mol );

};

// create canonical SMILES for the molecules, and sort using greater<string> and return
// paired with the original pointers in all_mols.
void create_smiles( std::vector<OEChem::OEMolBase *> &all_mols ,
                    std::vector<std::pair<std::string,OEChem::OEMolBase *> > &smiles );

// ****************************************************************************
class TooManyOutMols {

public :

  explicit TooManyOutMols( OEChem::OEMolBase &in_mol ) :
    in_mol_( in_mol ) {}

  OEChem::OEMolBase &in_mol_;

};

#endif // TAUTENUM_H
