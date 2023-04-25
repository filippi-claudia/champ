# PLT@NLeSC(2020)
# Ravindra Shinde (2021)
# (c) TREX-CoE (https://trex-coe.eu)

# Note:: How to use this ctest:
# 1. after a successful compilation, "cd build" and run "ctest"
# 2. To add more tests, just create a folder under CI_tests and have
#    a local bash scripts in that folder.
# 3. The following cmake ctest will automatically run all the tests
# 4. To just list the tests, run "ctest -N"

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()


file(GLOB_RECURSE tests_script LIST_DIRECTORIES false "*.sh")
  foreach(test_src ${tests_script})
      get_filename_component(dir ${test_src} PATH)
      get_filename_component(scriptname ${test_src} NAME)

      # Get filename and dir components as tags for the test
      string(REPLACE ".sh" "" fname_no_sh "${scriptname}")
      string(REPLACE "_" ";" tags "${fname_no_sh}")
      string(REPLACE "/" ";" dirs "${dir}")
      list(GET dirs -1 dir_tag)
      list(APPEND tags "${dir_tag}")

      # print the tags
      message("${BoldBlue} Adding test ${test_src} ${ColourReset}")
      message("   ${Blue}${tags}${ColourReset}")


      add_test(NAME ${test_src}
               WORKING_DIRECTORY ${dir}
               COMMAND "${scriptname}"
              )

      set_tests_properties("${test_src}" PROPERTIES LABELS "${tags}")

  endforeach(test_src)
