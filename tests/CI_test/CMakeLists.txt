# PLT@NLeSC(2020)

#message( "Binary:" ${CMAKE_BINARY_DIR} )
#message( "Source:" ${CMAKE_SOURCE_DIR} )
#message( "Current:" ${CMAKE_CURRENT_SOURCE_DIR} )
#message( "compiler:i" ${CMAKE_Fortran_COMPILER} )
#           ${CTEST_CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/${test_src})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CTEST_CMAKE_COMMAND ${CMAKE_BINARY_DIR}/vmc.mov1)

#  add_test(NAME free.inp 
#           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#           COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/script.sh 
#           ${CTEST_CMAKE_COMMAND} ${CMAKE_CURRENT_SOURCE_DIR}/free.inp)

file(GLOB tests_champ "*.inp")
  foreach(test_src ${tests_champ})
    add_test(NAME ${test_src}  
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/script.sh 
             ${CTEST_CMAKE_COMMAND} ${test_src})
  endforeach(test_src)

