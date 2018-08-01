# build the Utils library

add_library(Utils INTERFACE)
target_include_directories(Utils INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:include/utils>
 ) 

install(FILES complexify.h DESTINATION include/utils)
