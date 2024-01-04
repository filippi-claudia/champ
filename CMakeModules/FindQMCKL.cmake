#===========================================

# This file is distirbuted under the BSD 3-Clause License.
# Copyright (c) 2021, TREX Center of Excellence

#===========================================

message(" ")
message("Looking for the QMCKL library:")

set(QMCKL_SEARCH_PATHS
	~/Library/Frameworks
	/Library/Frameworks
	/usr/local
	/usr
	/sw # Fink
	/opt/local # DarwinPorts
	/opt/csw # Blastwave
	/opt
)

find_path(QMCKL_INCLUDE_DIR
	  NAMES qmckl_gpu.h
	  HINTS $ENV{QMCKL_DIR}
	  PATH_SUFFIXES include/qmckl include
	  PATHS ${QMCKL_SEARCH_PATHS}
	  )


# No need to specify platform-specific prefix (e.g. libqmckl on Unix) or
# suffix (e.g. .so on Unix or .dylib on MacOS) in NAMES. CMake takes care of that.
find_library(QMCKL_LIBRARY
             NAMES qmckl
	     HINTS $ENV{QMCKL_DIR}
	     PATH_SUFFIXES lib64 lib
	     PATHS ${QMCKL_SEARCH_PATHS}
	     )


# Handle the QUIETLY and REQUIRED arguments and set QMCKL_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(QMCKL DEFAULT_MSG QMCKL_LIBRARY QMCKL_INCLUDE_DIR )
MARK_AS_ADVANCED(QMCKL_INCLUDE_DIR QMCKL_LIBRARY)

# Mot setting _INCLUDE_DIR and _LIBRARIES is considered a bug,
# see https://gitlab.kitware.com/cmake/community/-/wikis/doc/tutorials/How-To-Find-Libraries
set(QMCKL_LIBRARIES ${QMCKL_LIBRARY})
set(QMCKL_INCLUDE_DIRS ${QMCKL_INCLUDE_DIR})
