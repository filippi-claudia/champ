function(gitinfo source)

  # Check is source is a valid path
  if(NOT EXISTS ${source})
    message(FATAL_ERROR "'${source}' is not a valid path")
  endif()

  # Define the possible location of the .git directory
  set(GIT_DIR "${source}/.git")

  # Check if .git folder exist
  if(EXISTS ${GIT_DIR})

    #
    set(GIT_DIR "${GIT_DIR}" CACHE PATH "Project .git directory")

    # Check if git is installed
    if(NOT GIT_FOUND)
      find_package(Git QUIET)
    endif()
    if(NOT GIT_FOUND)
      message(AUTHOR_WARNING "Git not found, cannot get git informations")
      return()
    endif()

    # whether or not the working tree is dirty
    execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain
            WORKING_DIRECTORY ${source}
            RESULT_VARIABLE exit_code
            OUTPUT_VARIABLE GIT_IS_DIRTY OUTPUT_STRIP_TRAILING_WHITESPACE)
    # the working tree is dirty when the error code is different from 0
    # or if the output is not empty
    if(NOT exit_code EQUAL 0 OR NOT ${GIT_IS_DIRTY} STREQUAL "")
      unset(GIT_IS_DIRTY)
      set(GIT_IS_DIRTY ON CACHE BOOL "Indicate if current branch is dirty")
    else()
      unset(GIT_IS_DIRTY)
      set(GIT_IS_DIRTY OFF CACHE BOOL "Indicate if current branch is dirty")
    endif()

    # name of the brack associated te HEAD
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_HEAD_BRANCH OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_HEAD_BRANCH "${GIT_HEAD_BRANCH}" CACHE INTERNAL "name of the brack associated te HEAD")

    # git revision full hash
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%H" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_REVISION_HASH OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_REVISION_HASH "${GIT_REVISION_HASH}" CACHE INTERNAL "git revision full hash")

    # shorten version of git revision
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%h" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_REVISION OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_REVISION "${GIT_REVISION}" CACHE INTERNAL "shorten version of git revision")

    # shorten version of git revision name
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%s" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_REVISION_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_REVISION_NAME "${GIT_REVISION_NAME}" CACHE INTERNAL "shorten version of git revision name")

    # author name
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%an" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_AUTHOR_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_AUTHOR_NAME "${GIT_AUTHOR_NAME}" CACHE INTERNAL "git author name")

    # author email
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%ae" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_AUTHOR_EMAIL OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_AUTHOR_EMAIL "${GIT_AUTHOR_EMAIL}" CACHE INTERNAL "git author email")

    # author date
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%ad" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_AUTHOR_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_AUTHOR_DATE "${GIT_AUTHOR_DATE}" CACHE INTERNAL "git author date")

    # author date, strict ISO 8601 format
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%aI" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_AUTHOR_DATE_ISO OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_AUTHOR_DATE_ISO "${GIT_AUTHOR_DATE_ISO}" CACHE INTERNAL "git author date, strict ISO 8601 format")

    # committer name
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%cn" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_COMMITTER_NAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_COMMITTER_NAME "${GIT_COMMITTER_NAME}" CACHE INTERNAL "git committer name")

    # committer email
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%ce" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_COMMITTER_EMAIL OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_COMMITTER_EMAIL "${GIT_COMMITTER_EMAIL}" CACHE INTERNAL "git committer email")

    # committer date
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%cd" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_COMMITTER_DATE OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_COMMITTER_DATE "${GIT_COMMITTER_DATE}" CACHE INTERNAL "git committer date")

    # committer date, strict ISO 8601 format
    execute_process(COMMAND ${GIT_EXECUTABLE} show -s "--format=%cI" HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_COMMITTER_DATE_ISO OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_COMMITTER_DATE_ISO "${GIT_COMMITTER_DATE_ISO}" CACHE INTERNAL "git committer date, strict ISO 8601 format")

    # origin remote url
    execute_process(COMMAND ${GIT_EXECUTABLE} config --get remote.origin.url
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_REMOTE_ORIGIN_URL OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_REMOTE_ORIGIN_URL "${GIT_REMOTE_ORIGIN_URL}" CACHE INTERNAL "git origin remote url")

    # most recent tag of the current branch
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --abbrev=0 HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_LATEST_TAG_LONG OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_LATEST_TAG_LONG "${GIT_LATEST_TAG_LONG}" CACHE INTERNAL "git most recent tag of the current branch")

    # most recent tagname of the current branch
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags HEAD
            WORKING_DIRECTORY ${source}
            OUTPUT_VARIABLE GIT_LATEST_TAG OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(GIT_LATEST_TAG "${GIT_LATEST_TAG}" CACHE INTERNAL "git most recent tagname of the current branch")

  ## Printing
  message("")
  message("")
  message("---------------------------------------------------------------------")
  message("Printing the version control information :: git")
  message("---------------------------------------------------------------------")
  message(STATUS "git dir                          :: ${GIT_DIR}")
  message("")
  message(STATUS "git head branch                  :: ${GIT_HEAD_BRANCH}")
  message(STATUS "git revision                     :: ${GIT_REVISION}")
  message(STATUS "git revision hash                :: ${GIT_REVISION_HASH}")
  message(STATUS "git revision name                :: ${GIT_REVISION_NAME}")
  message("")
  message(STATUS "git author name                  :: ${GIT_AUTHOR_NAME}")
  message(STATUS "git author email                 :: ${GIT_AUTHOR_EMAIL}")
  message(STATUS "git author date                  :: ${GIT_AUTHOR_DATE}")
  message(STATUS "git author date iso              :: ${GIT_AUTHOR_DATE_ISO}")
  message("")
  message(STATUS "git committer name               :: ${GIT_COMMITTER_NAME}")
  message(STATUS "git committer email              :: ${GIT_COMMITTER_EMAIL}")
  message(STATUS "git committer date               :: ${GIT_COMMITTER_DATE}")
  message(STATUS "git committer date iso           :: ${GIT_COMMITTER_DATE_ISO}")
  message("")
  message(STATUS "git remote origin url            :: ${GIT_REMOTE_ORIGIN_URL}")
  message(STATUS "git latest tag long              :: ${GIT_LATEST_TAG_LONG}")
  message(STATUS "git latest tag                   :: ${GIT_LATEST_TAG}")
  message("---------------------------------------------------------------------")
  else()
    message(STATUS "CHAMP source directory does not contain git information")
  endif()

endfunction()
