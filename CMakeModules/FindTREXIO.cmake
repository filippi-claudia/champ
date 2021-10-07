# If already in cache, be silent
if (TREXIO_INCLUDE_DIR AND TREXIO_LIBRARY)
  set (TREXIO_FIND_QUIETLY TRUE)
endif()

set(TREXIO_LIB "libtrexio.so")

find_path(TREXIO_INCLUDE_DIR NAMES trexio.h HINTS /usr/local/include)

find_library(TREXIO_LIBRARY
             NAMES ${TREXIO_LIB}
             PATHS /usr/lib
                   /usr/local/lib
                   ${HOME}/lib
             NO_DEFAULT_PATH)


# Handle the QUIETLY and REQUIRED arguments and set TREXIO_FOUND to TRUE if
# all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(TREXIO DEFAULT_MSG TREXIO_LIBRARY TREXIO_INCLUDE_DIR )
MARK_AS_ADVANCED(TREXIO_INCLUDE_DIR TREXIO_LIBRARY )
