file(GLOB PSPLINE_SOURCES LIST_DIRECTORIES false  *.f)

add_library(pspline ${PSPLINE_SOURCES})
set(PSPLINE_FLAGS ${Fortran_FLAGS})

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
list(REMOVE_ITEM PSPLINE_FLAGS "-fimplicit-none" "-Werror=implicit-interface")
endif()
target_compile_options(pspline PRIVATE ${PSPLINE_FLAGS})

