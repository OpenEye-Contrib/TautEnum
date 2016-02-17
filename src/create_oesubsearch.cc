//
// file create_oesubsearch.cc
// David Cosgrove
// AstraZeneca
// 26th November 2007
//
// This function takes a SMARTS string and returns an OESubSearch created on
// the stack from it. Throws a DACLIB::SMARTSDefnError exception if it can't,
// captured from the OpenEye function.

#include <string>
#include <vector>
#include <oechem.h>
#include "SMARTSExceptions.H"

using namespace std;
using namespace OEChem;
using namespace OEPlatform;
using namespace OESystem;

// **************************************************************************
namespace DACLIB {

  // **************************************************************************
  // reorder doesn't do anything at the moment but who knows, it might
  // one day.
  OESubSearch *create_oesubsearch( const string &smarts , bool reorder ) {

    oeosstream oeerrs;
    OEThrow.SetOutputStream( oeerrs );
    OESubSearch *subs = new OESubSearch( smarts.c_str() , reorder );
    OEThrow.SetOutputStream( oeerr );
    string errstr = oeerrs.str();
    if( !errstr.empty() ) {
      delete subs;
      subs = 0;
      // re-connect the OEThrow output stream
      OEThrow.SetOutputStream( OEPlatform::oeerr );
      throw SMARTSDefnError( errstr.c_str() );
    }

    // re-connect the OEThrow output stream
    OEThrow.SetOutputStream( OEPlatform::oeerr );
    return subs;

  }

  // **************************************************************************
  // this one uses sub_defns to expand smarts then calls the first one
  // sub_defns are pairs of name and SMARTS in that order.
  OESubSearch *create_eosubsearch( const string &smarts ,
                                   const string &smarts_name ,
                                   bool reorder ,
                                   vector<pair<string,string> > &sub_defns ) {

    string exp_smarts( smarts );
    if( string::npos != exp_smarts.find( "$") &&
        !OESmartsLexReplace( exp_smarts , sub_defns ) ) {
      exp_smarts.clear();
      throw( SMARTSSubDefnError( smarts , smarts_name ) );
    }
    cout << smarts << " expanded to " << exp_smarts << endl;
    return create_oesubsearch( exp_smarts , reorder );

  }

} // end of namespace DACLIB
