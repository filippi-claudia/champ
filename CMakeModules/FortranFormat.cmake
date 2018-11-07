
# Add fix or free format compiler flag
function (formatFortran SOURCE_FILES FLAG)

    # Intel does not care about the fixed format
  if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    foreach(_source ${SOURCE_FILES})
      set_source_files_properties(${_source} PROPERTIES FORTRAN_FORMAT ${FLAG})
      if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
	message (STATUS "${_source} is ${FLAG} format")
      endif (CMAKE_BUILD_TYPE EQUAL "DEBUG")

    endforeach()
  endif()
set(${SOURCE_FILES} ${SOURCE_FILES} PARENT_SCOPE)
endfunction()
