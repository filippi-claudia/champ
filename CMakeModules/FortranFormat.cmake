
# Add fix or free format compiler flag
function (formatFortran SOURCE_FILES FLAG)
  message(STATUS "Setting source format properties")

    # Intel does not care about the fixed format
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    foreach(_source ${SOURCE_FILES})
      set_source_files_properties(${_source} PROPERTIES FORTRAN_FORMAT ${FLAG})
      message (STATUS "${_source} is ${FLAG} format")
    endforeach()
  endif()
set(${SOURCE_FILES} ${SOURCE_FILES} PARENT_SCOPE)
endfunction()
