cmake_minimum_required(VERSION 3.4.0)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(ubitrack_core VERSION 1.3.0 LANGUAGES CXX)

include(${CMAKE_CURRENT_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

include(UbitrackConfig)

option(ENABLE_UNITTESTS "Enable Unit Tests" FALSE)
option(ENABLE_TRACING "Enable Tracing Providers" FALSE)

add_subdirectory(src)

if(ENABLE_UNITTESTS)
  add_subdirectory(tests)
endif(ENABLE_UNITTESTS)

ubitrack_install_doc("doc/bibliography.bib")













