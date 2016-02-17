cmake_minimum_required(VERSION 2.8)

project(TAUT_ENUM)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build_utils/CMakeModules/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -W -Wall  -Wshadow -Wno-error=shadow -Wno-error=unused-parameter -Wno-error=unused-function -fomit-frame-pointer -ffast-math -fvisibility=hidden -fvisibility-inlines-hidden -fPIC  -Werror " )

if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^3.1" )
  message( "Hopefully, the Ubuntu 14.04 machine, with a modern compiler." )
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
endif()
if( ${CMAKE_HOST_SYSTEM_VERSION} MATCHES "^2.6." )
  message( "It looks like a trusty old Centos 6 machine with an ancient compiler." )
endif()

find_package(OEToolkits COMPONENTS oedepict oechem oesystem oeplatform)
find_package(Boost COMPONENTS program_options regex iostreams filesystem system date_time thread REQUIRED)
find_package(Qt5 COMPONENTS Core Widgets REQUIRED)
find_package(Cairo REQUIRED)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

include_directories(${TAUT_ENUM_SOURCE_DIR}/inc
  ${TAUT_ENUM_SOURCE_DIR}/DACLib/inc)
include_directories( SYSTEM
  ${OEToolkits_INCLUDE_DIRS}
  ${Qt5Widgets_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS})

set(EXECUTABLE_OUTPUT_PATH ${TAUT_ENUM_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH ${TAUT_ENUM_SOURCE_DIR}/exe_${CMAKE_BUILD_TYPE})

# general library functions, that can be built into other programs, python modules etc.
set(TAUT_ENUM_LIB_SRCS
TautEnum.cc
TautStand.cc
smirks_helper_fns.cc
canned_tautenum_routines.cc)

# taut_enum - program for enumerating tautomers
set(TAUT_ENUM_SRCS
taut_enum.cc
TautEnumCallableBase.cc
TautEnumSettings.cc)

set(TAUT_ENUM_INCS
TautEnum.H
TautEnumCallableBase.H
TautEnumCallableSerial.H
TautEnumCallableThreaded.H
TautEnumSettings.H
TautStand.H
taut_enum_default_vector_bindings.H
taut_enum_default_enum_smirks_orig.H
taut_enum_default_enum_smirks_extended.H
taut_enum_default_standardise_smirks.H
taut_enum_protonate_a.H
taut_enum_protonate_b.H
taut_enum_protonate_vb.H)

set(TAUT_ENUM_DACLIB_SRCS
apply_daylight_arom_model_to_oemol.cc
build_time.cc
check_oechem_licence.cc
create_cansmi.cc
create_oesubsearch.cc
extract_smarts_from_smirks.cc
radical_atoms.cc)

set(TAUT_ENUM_DACLIB_INCS
chrono.h
SMARTSExceptions.H
)

# library for other building in to other stuff
add_library(tautenum STATIC ${TAUT_ENUM_LIB_SRCS} ${TAUT_ENUM_DACLIB_SRCS})

set(TAUT_ENUM_LIBS ${LIBS}
  ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES})

add_executable(taut_enum ${TAUT_ENUM_SRCS}
  ${TAUT_ENUM_INCS} ${TAUT_ENUM_DACLIB_SRCS} ${TAUT_ENUM_DACLIB_INCS})
target_link_libraries(taut_enum z tautenum ${TAUT_ENUM_LIBS} z pthread rt)

# mol_diff_viewer, for comparing files of SMILES string, e.g. from results
# of different tautomer enumerations

set(MOL_DIFF_VIEWER_SRCS
mol_diff_viewer.cc
MolDiffViewer.cc)

set(MOL_DIFF_VIEWER_INCS MolDiffViewer.H)

set(MOL_DIFF_VIEWER_DACLIB_SRCS
build_time.cc
check_oechem_licence.cc
create_cansmi.cc
draw_oemol_to_qimage.cc
split_smiles_into_atom_bits.cc
QTMolDisplay2D.cc
create_oesubsearch.cc
extract_smarts_from_smirks.cc)

set(MOL_DIFF_VIEWER_DACLIB_INCS
stddefs.H
DACOEMolAtomIndex.H
DACOEMolBondIndex.H
FileExceptions.H
QTMolDisplay2D.H
boost_tuples_and_bind.H)

set(MOL_DIFF_VIEWER_LIBS ${LIBS}
  ${Qt5Widgets_LIBRARIES}
  ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES} ${CAIRO_LIBRARY})

add_executable( mol_diff_viewer ${MOL_DIFF_VIEWER_SRCS} ${MOL_DIFF_VIEWER_DACLIB_SRCS}
  ${MOL_DIFF_VIEWER_INCS} ${MOL_DIFF_VIEWER_DACLIB_INCS} )
target_link_libraries(mol_diff_viewer ${MOL_DIFF_VIEWER_LIBS} z pthread rt)

# mol_diff_viewer2, for comparing files of SMILES string, e.g. from results
# of different tautomer enumerations. A 2nd version

set(MOL_DIFF_VIEWER2_SRCS
mol_diff_viewer2.cc
MolDiffViewer2.cc)

set(MOL_DIFF_VIEWER2_INCS MolDiffViewer2.H)

set(MOL_DIFF_VIEWER2_DACLIB_SRCS
build_time.cc
check_oechem_licence.cc
create_cansmi.cc
draw_oemol_to_qimage.cc
split_smiles_into_atom_bits.cc
QTMolDisplay2D.cc
create_oesubsearch.cc
extract_smarts_from_smirks.cc)

set(MOL_DIFF_VIEWER2_DACLIB_INCS
stddefs.H
DACOEMolAtomIndex.H
DACOEMolBondIndex.H
FileExceptions.H
QTMolDisplay2D.H
boost_tuples_and_bind.H)

set(MOL_DIFF_VIEWER2_LIBS ${LIBS}
  ${Qt5Widgets_LIBRARIES}
  ${OEToolkits_LIBRARIES}
  ${Boost_LIBRARIES} ${CAIRO_LIBRARY})

add_executable( mol_diff_viewer2 ${MOL_DIFF_VIEWER2_SRCS} ${MOL_DIFF_VIEWER2_DACLIB_SRCS}
  ${MOL_DIFF_VIEWER2_INCS} ${MOL_DIFF_VIEWER2_DACLIB_INCS} )
target_link_libraries(mol_diff_viewer2 ${MOL_DIFF_VIEWER2_LIBS} z pthread rt)
