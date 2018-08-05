# CoDiPack doesn't come with a CMake project file, so create a CMakeLists.txt
# for it
# IMPORTED targets don't seem to work correctly, so use INTERFACE instead


# options
option(ENABLE_ADEPT "enable Adept" OFF)
option(ADEPT_DIR "path to Adept")

#TODO: add error message of ADEPT_DIR is not set

message("ADEPT_DIR = ${ADEPT_DIR}")

if(ENABLE_ADEPT)

  # see: https://stackoverflow.com/questions/8774593/cmake-link-to-external-library/41909627#41909627
  #find_library(Adept_lib adept PATHS ${ADEPT_DIR}/lib)
  #  find_path(Adept_header adept.h PATHS ${ADEPT_DIR}/include)

  #  link_directories(${ADEPT_DIR})
  #add_library(Adept INTERFACE)
  #target_link_libraries(Adept INTERFACE ${ADEPT_DIR}/lib/libadept)
  #target_include_directories(Adept INTERFACE ${ADEPT_DIR}/include)


  #add_library(Adept SHARED IMPORTED)
  #set_target_properties(Adept PROPERTIES
  #  IMPORTED_LOCATION ${ADEPT_DIR}/lib/libadept
  #  INTERFACE_INCLUDE_DIRECTORIES ${ADEPT_DIR}/include
  #  INTERFACE_LINK_LIBRARIES ${ADEPT_DIR}/lib
  #)
endif()


