//
// file check_oechem_licence.cc
// David Cosgrove
// AstraZeneca
// 18th December 2007
//
// Checks the presence of an OEChem licence file pointed to by OE_LICENSE,
// makes an error message, and returns true or false.
// Uses the boost filesystem library to translate the file pointed to by
// OE_LICENSE and see if it exists, and the boost date_time library to
// translate the expiry date, if any, into a more human readable form.
// Bit of an overkill, that last bit, probably.

#include <iostream>
#include <string>
#include <boost/filesystem.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <oechem.h>

using namespace std;

// ****************************************************************************
namespace DACLIB {

  bool check_oechem_licence( string &err_msg ) {

    unsigned int expdate[3];

    if( !getenv( "OE_LICENSE" ) ) {
      err_msg =
          string( "No OpenEye licence defined.\nEnvironment variable OE_LICENSE\n\
                  must point to a valid licence file.\nCan't continue." );
    } else {
      if( !boost::filesystem::exists( getenv( "OE_LICENSE" ) ) ) {
        err_msg =
            string( "OpenEye Licence file " ) + string( getenv( "OE_LICENSE" ) ) +
            string( "\ndoesn't exist. Can\'t continue." );
      } else {
        if( !OEChem::OEChemIsLicensed( 0 , expdate ) ) {
          if( 0 == expdate[0] ) {
            err_msg = string( "No OEChem Licence found.\nCan\'t continue." );
          } else {
            boost::gregorian::date d( expdate[2] , expdate[1] , expdate[0] );
            err_msg =
                string( "Your OEChem Licence expired on " ) +
                boost::gregorian::to_simple_string( d ) +
                string( "\nCan\'t continue." );
          }
        }
      }
    }

    return err_msg.empty();

  }

} // EO namespace DACLIB
