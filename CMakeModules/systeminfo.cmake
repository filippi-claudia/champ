foreach(key
  IN ITEMS
    OS_NAME
    OS_RELEASE
    OS_VERSION
    OS_PLATFORM
    FQDN
    PROCESSOR_DESCRIPTION
    NUMBER_OF_LOGICAL_CORES
    NUMBER_OF_PHYSICAL_CORES
    TOTAL_VIRTUAL_MEMORY
    AVAILABLE_VIRTUAL_MEMORY
    TOTAL_PHYSICAL_MEMORY
    AVAILABLE_PHYSICAL_MEMORY
  )
  cmake_host_system_information(RESULT _${key} QUERY ${key})
endforeach()

message("")
message("")
message("---------------------------------------------------------------------")
message(" System Information ")
message("---------------------------------------------------------------------")
message(STATUS   "OS Name                          :: " ${_OS_NAME})
message(STATUS   "OS Release                       :: " ${_OS_RELEASE})
message(STATUS   "OS Version                       :: " ${_OS_VERSION})
message(STATUS   "OS Platform                      :: " ${_OS_PLATFORM})
message(STATUS   "Domain Name                      :: " ${_FQDN})
message(STATUS   "Processor                        :: " ${_PROCESSOR_DESCRIPTION})
message(STATUS   "Number of Logical Cores          :: " ${_NUMBER_OF_LOGICAL_CORES})
message(STATUS   "Number of Physical Cores         :: " ${_NUMBER_OF_PHYSICAL_CORES})
message(STATUS   "Total Virtual Memory             :: " ${_TOTAL_VIRTUAL_MEMORY} " MB")
message(STATUS   "Available Virtual Memory         :: " ${_AVAILABLE_VIRTUAL_MEMORY} " MB")
message(STATUS   "Total Physical Memory            :: " ${_TOTAL_PHYSICAL_MEMORY} " MB")
message(STATUS   "Available Physical Memory        :: " ${_AVAILABLE_PHYSICAL_MEMORY} " MB")
message("")

