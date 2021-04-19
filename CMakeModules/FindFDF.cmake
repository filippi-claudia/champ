
# If already in cache, be silent
if (FDF_INCLUDE_DIRS AND FDF_LIBRARIES)
  set (FDF_FIND_QUIETLY TRUE)
endif()

if(NOT BUILD_SHARED_LIBS)
  set(FDF_LIB "libfdf.a")
  set(FDF_INCLUDE "/usr/local/include")  
endif()

find_path(FDF_INCLUDE_DIR 
             NAMES ${FDF_INCLUDE}
             HINTS /usr/local/include             
             PATHS /usr/local/include
             NO_DEFAULT_PATH)

find_library(FDF_LIBRARY
             NAMES ${FDF_LIB}
             PATHS /usr/local/lib
             HINTS /usr/local/lib
             NO_DEFAULT_PATH)   

set(FDF_INCLUDE_DIRS ${FDF_INCLUDE})
set(FDF_LIBRARIES ${FDF_LIBRARY})


INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FDF DEFAULT_MSG FDF_LIBRARIES FDF_INCLUDE_DIRS)

MARK_AS_ADVANCED(FDF_INCLUDE_DIRS FDF_LIBRARIES)
