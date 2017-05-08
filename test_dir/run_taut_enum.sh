#!/bin/bash

# standardise 1,2,4-triazoles
../src/exe_DEBUG/taut_enum -I triazoles.smi -O triazoles_std.smi --standardise-only

# standardise 1,2,4-triazoles and do original enumeration
../src/exe_DEBUG/taut_enum -I triazoles.smi -O triazoles_orig_enum.smi --original-enumeration

# standarise 1,2,4-triazoles and do extended enumeration
../src/exe_DEBUG/taut_enum -I triazoles.smi -O triazoles_ext_enum.smi --extended-enumeration

# standarise 1,2,4-triazoles and do both enumerations
../src/exe_DEBUG/taut_enum -I triazoles.smi -O .smi --original-enumeration | \
    ../src/exe_DEBUG/taut_enum -I .smi -O triazoles_full_enum.smi --extended-enumeration
