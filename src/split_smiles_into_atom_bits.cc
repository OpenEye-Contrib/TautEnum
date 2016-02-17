//
// file split_smiles_into_atom_bits.cc
// David Cosgrove
// AstraZeneca
// 18th February 2008
//
// Takes a SMILES string and splits it up into individual portions, one portion
// per atom. Doesn't include extra stuff like ( for branching, ring closure
// numbers etc. and also doesn't deal with 2 character atomic symbols for
// anything other than Br and Cl. I think that's ok for what it's going to be
// used for.

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace std;

// ****************************************************************************
namespace DACLIB {

// starting from position i, assumed to contain a [, move through
// smiles to the next ] that closes opening one, allowing for the fact
// that there may be nested [/] pairs. That's only if smiles is, in fact
// a recursive SMARTS rather than a SMILES, but split_smiles_into_atom_bits
// is used for that, too.
string extract_to_close_square_bracket( const string &smiles ,
                                        unsigned int &i ) {

  string retval;
  retval += smiles[i++];
  int brac_depth = 1;
  while( i < smiles.length() ) {
    retval += smiles[i];
    switch( smiles[i] ) {
    case '[' :
      ++brac_depth;
      break;
    case ']' :
      --brac_depth;
      if( !brac_depth ) {
        return retval;
      }
      break;
    }
    ++i;
  }

  // if we're here, there's probably an error with the SMILES, but we won't
  // worry about that.
  return retval;

}

// ****************************************************************************
void split_smiles_into_atom_bits( const string &smiles ,
                                  vector<string> &bits ) {

  unsigned int i = 0;
  while( i < smiles.length() ) {
    bits.push_back( string( "" ) );
    switch( smiles[i] ) {
    case '[' :
      // extract... increments i as it goes along
      bits.back() = extract_to_close_square_bracket( smiles , i );
      break;
    case 'C' :
      bits.back() += 'C';
      if( 'l' == smiles[i+1] ) {
        bits.back() += 'l';
        ++i;
      }
      ++i;
      break;
    case 'B' :
      bits.back() += 'B';
      if( 'r' == smiles[i+1] ) {
        bits.back() += 'r';
        ++i;
      }
      ++i;
      break;
    default :
      bits.back() += smiles[i];
      ++i;
      break;
    }
    while( i < smiles.length() && '[' != smiles[i] && '*' != smiles[i] &&
           !isalpha( smiles[i] ) ) {
      ++i;
    }
  }

#ifdef NOTYET
  cout << "split : " << smiles << endl;
  copy( bits.begin() , bits.end() ,
        ostream_iterator<string>( cout , " " ) );
  cout << endl;
#endif

}

} // end of namespace

