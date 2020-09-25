# PLT@NLeSC(2020)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(CTEST_CMAKE_COMMAND ${CMAKE_BINARY_DIR}/vmc.mov1)

file(GLOB tests_champ "*.inp")
  foreach(test_src ${tests_champ})
    add_test(NAME ${test_src}  
             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
             COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/script.sh 
             ${CTEST_CMAKE_COMMAND} ${test_src})
  endforeach(test_src)

