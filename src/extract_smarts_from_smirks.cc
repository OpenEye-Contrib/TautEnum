//
// file extract_smarts_from_smirks.cc
// David Cosgrove
// AstraZeneca
// 5th December 2011
//
// Crude, text-based extraction of a SMARTS string from a SMIRKS.

#include <string>

using namespace std;

namespace DACLIB {

// *******************************************************************************
string extract_smarts_from_smirks( const string &smirks_string ) {

  // first, extract the SMARTS from the SMIRKS and tidy up any explicit H atoms
  // which cannot be used
  // split the first half of the SMIRKS out
  size_t i = smirks_string.find( ">>" );
  string smarts_string = smirks_string.substr( 0 , i );
  if( string::npos != i ) {
    // the SMIRKS could well have [H] strings in them, which can't be in a SMARTS
    // so take them out, including [H:3] type ones as well
    i = smarts_string.find( "[H" );
#ifdef NOTYET
    cout << "smarts string " << smarts_string << endl;
#endif
    while( i != string::npos ) {
      size_t j = smarts_string.find( ']' , i );
      for( size_t k = i ; k <= j ; ++k ) {
        smarts_string[k] = ' ';
      }
      smarts_string.replace( i , j - i + 1 , "" );
      i = smarts_string.find( "[H" );
    }
    // now remove any () that might be left
#ifdef NOTYET
    cout << "smarts string " << smarts_string << endl;
#endif
    i = smarts_string.find( "()" );
    while( i != string::npos ) {
      smarts_string.replace( i , 2 , "" );
      i = smarts_string.find( "()" );
    }
#ifdef NOTYET
    cout << "smarts string " << smarts_string << endl;
#endif
  }

  return smarts_string;

}

} // EO namespace DACLIB
