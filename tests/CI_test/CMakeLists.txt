# PLT@NLeSC(2020)
# Ravindra Shinde (2021)
# (c) TREX-CoE (https://trex-coe.eu)

# Note:: How to use this ctest:
# 1. after a successful compilation, "cd build" and run "ctest"
# 2. To add more tests, just create a folder under CI_tests and have
#    a local bash scripts in that folder.
# 3. The following cmake ctest will automatically run all the tests
# 4. To just list the tests, run "ctest -N"

file(GLOB_RECURSE tests_script LIST_DIRECTORIES false "*.sh")
  foreach(test_src ${tests_script})
      get_filename_component(dir ${test_src} PATH)
      get_filename_component(scriptname ${test_src} NAME)
      add_test(NAME ${test_src}
               WORKING_DIRECTORY ${dir}
              COMMAND ${scriptname}
              )
  endforeach(test_src)

