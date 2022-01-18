# PLT@NLeSC(2020)
# Ravindra Shinde (2021)
# (c) TREX-CoE (https://trex-coe.eu)

# Note:: How to use this ctest:
# 1. after a successful compilation, "cd build" and run "ctest"
# 2. To add more tests, just create a folder under CI_tests and have
#    a local bash scripts in that folder.
# 3. The following cmake ctest will automatically run all the tests
# 4. To just list the tests, run "ctest -N"

file(GLOB_RECURSE tests_script LIST_DIRECTORIES false "*_short.sh")
  foreach(test_src ${tests_script})
      get_filename_component(dir ${test_src} PATH)
      get_filename_component(scriptname ${test_src} NAME)
      string(REPLACE ".sh" "" tags "${scriptname}")
      string(REPLACE "_" ";" tags "${tags}")
      string(REPLACE "/" ";" dirs "${dir}")
      list(GET dirs -1 dir_tag)
      list(APPEND tags "${dir_tag}")

      add_test(NAME "${test_src}.1"
               WORKING_DIRECTORY ${dir}
               COMMAND "bash" "-c" 
               "./${scriptname} > ${scriptname}.1.out 2> >(tee ${scriptname}.1.err >&2)"
              )
      set_tests_properties("${test_src}.1" PROPERTIES LABELS "${tags};one_proc")
      add_test(NAME "${test_src}.2"
               WORKING_DIRECTORY ${dir}
               COMMAND "bash" "-c" 
               "mpirun -np 2 ./${scriptname} > ${scriptname}.2.out 2> >(tee ${scriptname}.2.err >&2)"
              )
      set_tests_properties("${test_src}.2" PROPERTIES LABELS "${tags};two_proc")
  endforeach(test_src)

