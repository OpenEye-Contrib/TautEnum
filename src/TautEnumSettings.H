//
// file TautEnumSettings.H
// David Cosgrove
// AstraZeneca
// 31st August 2011
//
// Settings interface for program taut_enum

#ifndef TAUTENUMSETTINGS_H
#define TAUTENUMSETTINGS_H

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

// **************************************************************************

class TautEnumSettings {

public :

  TautEnumSettings( int argc , char **argv );

  void print_usage( std::ostream &os ) const;
  void print_error( std::ostream &os ) const;

  std::string input_mol_file() const { return in_mol_file_; }
  std::string output_mol_file() const { return out_mol_file_; }
  std::string vb_file() const { return vb_file_; }
  std::string standardise_smirks_file() const { return stand_smirks_file_; }
  std::string enumerate_smirks_file() const { return enum_smirks_file_; }
  std::string protonation_standardisation_file() const { return prot_stand_file_; }
  std::string protonation_enumeration_file() const { return prot_enum_file_; }
  std::string protonation_vb_file() const { return prot_vb_file_; }
  std::string name_postfix() const { return name_postfix_; }
  bool standardise_only() const { return standardise_only_; }
  bool original_enumeration() const { return orig_enumeration_; }
  bool extended_enumeration() const { return extended_enumeration_; }
  bool add_numbers_to_name() const { return add_numbers_to_name_; }
  bool add_smirks_to_name() const { return add_smirks_to_name_; }
  bool canonical_tautomer() const { return canon_taut_; }
  bool enumerate_protonation() const { return enum_protonation_; }
  bool include_input_in_output() const { return inc_input_in_output_; }
  bool strip_salts() const { return strip_salts_; }
  unsigned int max_tautomers() const { return max_tauts_; }
  bool do_threaded() const { return do_threaded_; }
  int num_threads() const { return num_threads_; } // -1 means use all available threads
  bool verbose() const { return verbose_; }

  bool operator!() const;

private :

  std::string in_mol_file_ , out_mol_file_;
  std::string vb_file_;
  std::string enum_smirks_file_;
  std::string stand_smirks_file_;
  std::string prot_stand_file_;
  std::string prot_enum_file_;
  std::string prot_vb_file_;
  std::string name_postfix_; // postfix for name before number when add_numbers_to_name_ is true. Defaults to _
  mutable bool standardise_only_;
  bool orig_enumeration_;
  bool extended_enumeration_;
  bool add_numbers_to_name_;
  bool add_smirks_to_name_;
  bool canon_taut_;
  bool enum_protonation_;
  bool inc_input_in_output_;
  bool strip_salts_; // whether to keep just largest component or not
  unsigned int max_tauts_;
  bool do_threaded_;
  int num_threads_;
  bool verbose_;

  std::string usage_text_;
  mutable std::string error_msg_;

  void build_program_options( boost::program_options::options_description &desc );

};


#endif // TAUTENUMSETTINGS_H
