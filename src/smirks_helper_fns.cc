//
// file smirks_helper_fns.cc
// David Cosgrove
// AstraZeneca
// 2nd August 2011
//
// Some functions of use to various SMIRKS programs

#include "FileExceptions.H"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include <oechem.h>

using namespace std;
using namespace OEChem;

typedef boost::shared_ptr<OEChem::OELibraryGen> pOELibGen;
typedef boost::shared_ptr<OEChem::OEMolBase> pOEMol;

namespace DACLIB {

// ************************************************************************************
// OESmartsLexReplace messes up if there are duplicate vector bindings,
// even if they're for the same thing.  check_for_duplicate_vector_bindings()
// checks for this. It issues a warning if they're the same name and definition,
// and removes one of the instances, emits an error message and aborts if
// the definitions are different.
void check_for_duplicate_vector_bindings( vector<pair<string,string> > &vbs ) {

  using namespace boost;
  sort( vbs.begin() , vbs.end() ,
        bind( less<string>() ,
              bind( &pair<string,string>::first , _1 ) ,
              bind( &pair<string,string>::first , _2 ) ) );
  vector<pair<string,string> >::iterator p = vbs.begin();
  while( 1 ) {
    p = adjacent_find( vbs.begin() , vbs.end() ,
                       bind( equal_to<string>() ,
                             bind( &pair<string,string>::first , _1 ) ,
                             bind( &pair<string,string>::first , _2 ) ) );
    if( p == vbs.end() ) {
      break;
    }
    if( p->second == (p+1)->second ) {
      cerr << "Warning duplicate definition for SMARTS vector binding "
           << p->first << " defined as " << p->second << endl;
      vbs.erase( p );
    } else {
      cerr << "Error : duplicate names for different vector bindings : " << endl
           << p->first << " defined as " << p->second
           << " and " << (p+1)->second << endl << "Program aborts." << endl;
      exit( 1 );
    }
  }

}

// ************************************************************************************
void read_vbs_from_istream( istream &is ,
                            vector<pair<string,string> > &vbs ) {

  string next_line;
  while( 1 ) {
    next_line = "";
    getline( is , next_line );
    if( ( is.eof() || !is.good() ) && next_line.empty() ) {
      break;
    }
    if( '#' == next_line[0] ) {
      continue; // lines that start with # are comments
    }
    vector<string> splits;
    boost::algorithm::trim( next_line );
    if( next_line.empty() ) {
      continue; // trimmed down to nothing
    }
    boost::algorithm::split( splits , next_line , boost::algorithm::is_any_of( " \t" ) , boost::algorithm::token_compress_on );
    if( splits.size() < 2 ) {
      continue; // something's wrong, though
    }
    vbs.push_back( make_pair( splits[0] , splits[1] ) );
  }

  // cout << "Read " << vbs.size() << " vector bindings." << endl;

  // OESmartsLexReplace at the moment (version 1.7.2.4, 16th May 2011) gives rubbish
  // back if a vector binding is defined more than once, so check for that.
  check_for_duplicate_vector_bindings( vbs );

}

// ************************************************************************************
void read_vbs_from_file( const string &vb_file ,
                         vector<pair<string,string> > &vbs ) {

  if( !vb_file.empty() ) {
    ifstream ifs( vb_file.c_str() );
    if( !ifs || !ifs.good() ) {
      throw DACLIB::FileReadOpenError( vb_file.c_str() );
    }

    read_vbs_from_istream( ifs , vbs );
  }

}

// ************************************************************************************
void read_vbs_from_string( const string &vbs_string ,
                         vector<pair<string,string> > &vbs ) {

  istringstream iss( vbs_string );
  read_vbs_from_istream( iss , vbs );

}

// ****************************************************************************
void read_smirks_from_istream( istream &is , vector<pair<string,string> > &smks ) {

  string next_line;
  while( 1 ) {
    next_line = "";
    getline( is , next_line );
    if( ( is.eof() || !is.good() ) && next_line.empty() ) {
      break;
    }
    if( '#' == next_line[0] ) {
      continue; // lines that start with # are comments
    }
    vector<string> splits;
    boost::algorithm::trim( next_line );
    if( next_line.empty() ) {
      continue; // trimmed down to nothing
    }
    boost::algorithm::split( splits , next_line , boost::algorithm::is_any_of( " \t" ) , boost::algorithm::token_compress_on );
    if( 1 == splits.size() ) {
      smks.push_back( make_pair( string( "Smk" ) + boost::lexical_cast<string>( smks.size() + 1 ) ,
                                 splits[0] ) );
    } else {
      smks.push_back( make_pair( splits[1] , splits[0] ) );
    }
  }

  // cout << "Read " << smks.size() << " SMIRKS definitions." << endl;

 }

// ****************************************************************************
void read_smirks_from_file( const string &smirks_file ,
                            vector<pair<string,string> > &smks ) {

  ifstream ifs( smirks_file.c_str() );
  if( !ifs || !ifs.good() ) {
    throw DACLIB::FileReadOpenError( smirks_file.c_str() );
  }

  read_smirks_from_istream( ifs , smks );

}

// ****************************************************************************
void read_smirks_from_string( const string &smirks_string ,
                              vector<pair<string,string> > &smks ) {

  istringstream iss( smirks_string );

  read_smirks_from_istream( iss , smks );

}

// ************************************************************************************
void expand_vector_bindings( const vector<pair<string,string> > &in_smirks ,
                             vector<pair<string,string> > &vbs ,
                             vector<string> &exp_smirks ) {

  exp_smirks.clear();
  for( size_t i = 0 , is = in_smirks.size() ; i < is ; ++i ) {
    string tmp = in_smirks[i].second;
    if( !vbs.empty() ) {
      OESmartsLexReplace( tmp , vbs );
    }
    exp_smirks.push_back( tmp );
#ifdef NOTYET
    cout << "SMIRKS name = " << in_smirks[i].first << " value = " << in_smirks[i].second
         << " expands to " << tmp << endl;
#endif
  }

}

// ************************************************************************************
// make an OELibraryGen obect from the SMIRKS and set it up how we like it
pOELibGen create_libgen( const string &smirks ) {

  pOELibGen ret_val( new OEChem::OELibraryGen( smirks.c_str() ) );
  ret_val->SetAssignMapIdx( true );
  ret_val->SetExplicitHydrogens( true );

  return ret_val;

}

// ************************************************************************************
// create a set of OELibaryGen objects from the input SMIRKS
void create_libgens( const vector<string> &exp_smirks ,
                     const vector<pair<string,string> > &in_smirks ,
                     vector<pOELibGen> &lib_gens ) {

  for( size_t i = 0 , is = exp_smirks.size() ; i < is ; ++i ) {
    lib_gens.push_back( DACLIB::create_libgen( exp_smirks[i] ) );
    if( !*lib_gens.back() ) {
      cerr << "AWOOGA : error parsing SMIRKS " << exp_smirks[i]
              << " built from " << in_smirks[i].first << " : " << in_smirks[i].second << endl;
      exit( 1 );
    }
  }

}

} // EO namespace DACLIB
