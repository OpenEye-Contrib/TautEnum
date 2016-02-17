//
// file TautEnumSettings.cc
// David Cosgrove
// AstraZeneca
// 31st August 2011
//

#include "TautEnumSettings.H"

#include <iostream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

using namespace std;
namespace po = boost::program_options;

// ********************************************************************************
TautEnumSettings::TautEnumSettings( int argc , char **argv ) :
  name_postfix_( "_" ) , standardise_only_( true ) ,
  orig_enumeration_( false ) , extended_enumeration_( false ) ,
  add_numbers_to_name_( false ) , add_smirks_to_name_( false ) ,
  canon_taut_( false ) , enum_protonation_( false ) ,
  inc_input_in_output_( false ) , strip_salts_( false ) , max_tauts_( 256 ) ,
  do_threaded_( false ) , num_threads_( -1 ) , verbose_( false ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ********************************************************************************
void TautEnumSettings::print_usage( std::ostream &os ) const {

  os << usage_text_ << endl;

}

// ********************************************************************************
void TautEnumSettings::print_error( std::ostream &os ) const {

  os << error_msg_ << endl;

}

// ********************************************************************************
bool TautEnumSettings::operator!() const {

  if( orig_enumeration_ || extended_enumeration_ || enum_protonation_ || !enum_smirks_file_.empty() ) {
    standardise_only_ = false;
  }

  if( orig_enumeration_ && extended_enumeration_ ) {
    error_msg_ = "You can't have both original and extended enumerations.";
    return true;
  }
  if( include_input_in_output() && canonical_tautomer() ) {
    error_msg_ = "You can't have both canonical tautomer and include input in output.";
    return true;
  }
  if( in_mol_file_.empty() ) {
    error_msg_ = "You must specify an input molecule file.";
    return true;
  }
  if( out_mol_file_.empty() ) {
    error_msg_ = "You must specify an output molecule file.";
    return true;
  }

  return false;

}

// ********************************************************************************
void TautEnumSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
      ( "help" , "Produce this help text." )
      ( "input-molecule-file,I" , po::value<string>( &in_mol_file_ ) ,
        "Input molecule filename" )
      ( "output-molecule-file,O" , po::value<string>( &out_mol_file_ ) ,
        "Output molecule filename" )
      ( "standardise-smirks-file,S" , po::value<string>( &stand_smirks_file_ ) ,
        "File of SMIRKS transformations for standardisations." )
      ( "standardize-smirks-file" , po::value<string>( &stand_smirks_file_ ) ,
        "File of SMIRKS transformations for standardisations." )
      ( "enumerate-smirks-file,E" , po::value<string>( &enum_smirks_file_ ) ,
        "File of SMIRKS transformations for enumerations." )
      ( "vector-bindings-file,V" , po::value<string>( &vb_file_ ) ,
        "Name of file of vector bindings." )
      ( "name-postfix" , po::value<string>( &name_postfix_ ) ,
        "Postfix to molecule name before tautomer number. Defaults to \'_\'.")
      ( "standardise-only" , po::value<bool>( &standardise_only_ )->zero_tokens() ,
        "Just put each molecule in standard tautomer." )
      ( "standardize-only" , po::value<bool>( &standardise_only_ )->zero_tokens() ,
        "Just put each molecule in standard tautomer." )
      ( "original-enumeration" , po::value<bool>( &orig_enumeration_ )->zero_tokens() ,
        "Limited enumeration, akin to the original Leatherface." )
      ( "extended-enumeration" , po::value<bool>( &extended_enumeration_ )->zero_tokens() ,
        "Extended enumeration, with full keto-enol tautomerisation, for example.")
      ( "enumerate-protonation" , po::value<bool>( &enum_protonation_ )->zero_tokens() ,
        "Enumerate protonation states for molecules." )
      ( "add-numbers-to-name" , po::value<bool>( &add_numbers_to_name_ )->zero_tokens() ,
        "Add name postfix and tautomer number to name of each output molecule." )
      ( "add-smirks-to-name" , po::value<bool>( &add_smirks_to_name_ )->zero_tokens() ,
        "Add to the name a space-separated list of the SMIRKS names that were used to generate this tautomer from the input structure." )
      ( "canonical-tautomer" , po::value<bool>( &canon_taut_ )->zero_tokens() ,
        "Just output the canonical tautomer of each molecule." )
      ( "include-input-in-output" , po::value<bool>( &inc_input_in_output_ )->zero_tokens() ,
        "Before each set of tautomers is output, write the input molecule as it came in.")
      ( "strip-salts" , po::value<bool>( &strip_salts_)->zero_tokens() ,
        "Strip out all but largest component before processing." )
      ( "max-tautomers" , po::value<unsigned int>( &max_tauts_ ) ,
        "Maximum number of tautomers per molecule." )
      ( "do-threaded" , po::value<bool>( &do_threaded_ )->zero_tokens() ,
        "Whether to do a threaded run as opposed to default serial." )
      ( "num-threads" , po::value<int>( &num_threads_ ) ,
        "Number of threads to use. A number <= 0 means subtract from hardware thread maximum." )
      ( "verbose" , po::value<bool>( &verbose_ )->zero_tokens() ,
        "Extra output saying what's been going on." )
      ( "warm-feeling,W" , po::value<bool>( &verbose_ )->zero_tokens() ,
        "Extra output saying what's been going on." );

}
