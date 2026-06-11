# CHAMP CI integration tests
# Ravindra Shinde (2026)
#
# Each test folder carries a declarative `test.json` manifest; the
# manifests are validated, staged and executed by champ_test_runner.py.
# This file registers one ctest per case.  See README.md in this
# directory for usage and the manifest reference.

find_package(Python3 COMPONENTS Interpreter)
if(NOT Python3_Interpreter_FOUND)
  message(WARNING "Python3 interpreter not found -- CHAMP CI tests are disabled")
  return()
endif()

set(CHAMP_TEST_RUNNER "${CMAKE_CURRENT_SOURCE_DIR}/champ_test_runner.py")

# The MPI launcher comes from FindMPI (already required at the top level)
# but can be overridden, e.g. for srun or to add --oversubscribe.
set(CHAMP_TEST_MPIEXEC "${MPIEXEC_EXECUTABLE}" CACHE STRING
    "MPI launcher used to run the CI tests")
set(CHAMP_TEST_NUMPROC_FLAG "${MPIEXEC_NUMPROC_FLAG}" CACHE STRING
    "Flag of CHAMP_TEST_MPIEXEC that sets the number of ranks")
set(CHAMP_TEST_MPIEXEC_PREFLAGS "${MPIEXEC_PREFLAGS}" CACHE STRING
    "Extra launcher flags placed before the executable (e.g. --oversubscribe)")

# Tests run the binaries this build produces unless pointed elsewhere
# (e.g. at a different build mode; all builds emit into <source>/bin).
set(CHAMP_TEST_VMC_BINARY "$<TARGET_FILE:vmc.mov1>" CACHE STRING
    "vmc.mov1 binary exercised by the CI tests (default: this build's)")
set(CHAMP_TEST_DMC_BINARY "$<TARGET_FILE:dmc.mov1>" CACHE STRING
    "dmc.mov1 binary exercised by the CI tests (default: this build's)")

# A quick sanity test of the runner itself (no CHAMP binaries needed).
add_test(NAME runner-selftest
         COMMAND ${Python3_EXECUTABLE} ${CHAMP_TEST_RUNNER} selftest)
set_tests_properties(runner-selftest PROPERTIES LABELS "selftest" TIMEOUT 120)

# Enumerate and validate all manifests in one configure-time call.
file(GLOB _champ_manifests CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*/test.json")
list(SORT _champ_manifests)
if(NOT _champ_manifests)
  message(WARNING "No test.json manifests found under ${CMAKE_CURRENT_SOURCE_DIR}")
  return()
endif()

# Re-run the configure step when a manifest (or the runner) is edited, so
# label/timeout/case changes are picked up by the next build automatically.
set_property(DIRECTORY APPEND PROPERTY
             CMAKE_CONFIGURE_DEPENDS ${_champ_manifests} ${CHAMP_TEST_RUNNER})

execute_process(
  COMMAND ${Python3_EXECUTABLE} ${CHAMP_TEST_RUNNER} list ${_champ_manifests}
  OUTPUT_VARIABLE _champ_cases
  ERROR_VARIABLE _champ_notes
  RESULT_VARIABLE _champ_rc)
if(NOT _champ_rc EQUAL 0)
  message(FATAL_ERROR "CI test manifest validation failed:\n${_champ_notes}")
endif()

set(_n_added 0)
set(_n_skipped 0)
set(_skipped "")

string(REPLACE "\n" ";" _champ_case_lines "${_champ_cases}")
foreach(_line IN LISTS _champ_case_lines)
  if(_line STREQUAL "")
    continue()
  endif()
  string(REPLACE "\t" ";" _fields "${_line}")
  list(LENGTH _fields _nfields)
  if(NOT _nfields EQUAL 7)
    message(FATAL_ERROR "Malformed runner output line: ${_line}")
  endif()
  list(GET _fields 0 _manifest)
  list(GET _fields 1 _folder)
  list(GET _fields 2 _case)
  list(GET _fields 3 _nproc)
  list(GET _fields 4 _timeout)
  list(GET _fields 5 _requires)
  list(GET _fields 6 _labels)

  # Skip tests whose required features are not part of this build.
  set(_unmet "")
  if(NOT _requires STREQUAL "-")
    string(REPLACE "," ";" _reqlist "${_requires}")
    foreach(_req IN LISTS _reqlist)
      if(_req STREQUAL "trexio" AND NOT TREXIO_FOUND)
        set(_unmet "${_req}")
      elseif(_req STREQUAL "qmckl" AND NOT QMCKL_FOUND)
        set(_unmet "${_req}")
      endif()
    endforeach()
  endif()
  if(NOT _unmet STREQUAL "")
    math(EXPR _n_skipped "${_n_skipped}+1")
    list(APPEND _skipped "${_folder}/${_case} (needs ${_unmet})")
    continue()
  endif()

  string(REPLACE "," ";" _labellist "${_labels}")
  add_test(NAME "${_folder}/${_case}"
           COMMAND ${Python3_EXECUTABLE} ${CHAMP_TEST_RUNNER} run
                   --manifest "${_manifest}"
                   --case "${_case}"
                   --scratch "${CMAKE_CURRENT_BINARY_DIR}/scratch/${_folder}/${_case}"
                   --vmc "${CHAMP_TEST_VMC_BINARY}"
                   --dmc "${CHAMP_TEST_DMC_BINARY}"
                   --mpiexec "${CHAMP_TEST_MPIEXEC}"
                   "--numproc-flag=${CHAMP_TEST_NUMPROC_FLAG}"
                   "--mpiexec-preflags=${CHAMP_TEST_MPIEXEC_PREFLAGS}"
           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties("${_folder}/${_case}" PROPERTIES
    LABELS "${_labellist}"
    PROCESSORS "${_nproc}"
    TIMEOUT "${_timeout}")
  math(EXPR _n_added "${_n_added}+1")
endforeach()

# Disabled manifests are reported by the runner on stderr as "DISABLED <folder>".
set(_disabled "")
string(REPLACE "\n" ";" _champ_note_lines "${_champ_notes}")
foreach(_line IN LISTS _champ_note_lines)
  if(_line MATCHES "^DISABLED\t(.+)$")
    list(APPEND _disabled "${CMAKE_MATCH_1}")
  endif()
endforeach()

message(STATUS "CHAMP CI tests                   :: ${_n_added} cases registered")
if(_n_skipped GREATER 0)
  message(STATUS "CHAMP CI tests skipped           :: ${_n_skipped} (missing features in this build)")
  foreach(_entry IN LISTS _skipped)
    message(STATUS "   skipped :: ${_entry}")
  endforeach()
endif()
foreach(_entry IN LISTS _disabled)
  message(STATUS "   disabled in manifest :: ${_entry}")
endforeach()
