
# Add fixed format property
macro (formatFortran SOURCE_FILES)
  foreach(_source ${SOURCE_FILES})
    set_property(SOURCE ${_source} PROPERTY Fortran_FORMAT FIXED)
  endforeach(_source)
endmacro()
