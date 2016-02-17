# attempts to find Cairo v2.0
# it will define
# CAIRO_FOUND it if finds the library it needs
# CAIRO_INCLUDE_DIR
# CAIRO_LIBRARY

FIND_PATH(CAIRO_INCLUDE_DIRS cairo.h
  /usr/local/include
  /usr/include
  /opt/include )

if( NOT ${CAIRO_INCLUDE_DIR} )
  message("CAIRO_INCLUDE_DIR found as ${CAIRO_INCLUDE_DIR}")
else( )
  message("CAIRO_INCLUDE_DIR not found.")
endif ()

IF (WIN32)
  SET(LIBCAIRO "cairo.lib")
ELSE (WIN32)
  SET(LIBCAIRO "libcairo.so.2")
ENDIF (WIN32)

IF(CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64")
FIND_FILE(CAIRO_LIBRARY ${LIBCAIRO}
  /usr/local/lib64
  /usr/lib64)
ELSE()
FIND_FILE(CAIRO_LIBRARY ${LIBCAIRO}
  /usr/local/lib
  /usr/lib)
ENDIF()


message( "CAIRO_LIBRARY : ${CAIRO_LIBRARY}" )

if( ${CAIRO_LIBRARY} )
  set(CAIRO_FOUND True)
else()
  set(CAIRO_FOUND "" )
endif()