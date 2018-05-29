cmake_minimum_required(VERSION 3.4.0)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(ubitrack_core VERSION 1.3.0 LANGUAGES CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(TARGETS)

include(UbitrackConfig)

option(ENABLE_UNITTESTS "Enable Unit Tests" FALSE)
option(ENABLE_TRACING "Enable Tracing Providers" FALSE)

add_subdirectory(src)

if(ENABLE_UNITTESTS)
  add_subdirectory(tests)
endif(ENABLE_UNITTESTS)

file(GLOB _doc_files LIST_DIRECTORIES false "doc/bibliography.bib")
foreach(pfile ${_doc_files})
  string(REGEX REPLACE "^.*/doc/" "" pfile2 "${pfile}")
  GET_FILENAME_COMPONENT(fpath ${pfile2} PATH)
  IF(fpath)
    install(FILES ${pfile} DESTINATION "${UBITRACK_DOC_INSTALL_DIRECTORY}/${fpath}" COMPONENT doc)
  ELSE(fpath)
    install(FILES ${pfile} DESTINATION "${UBITRACK_DOC_INSTALL_DIRECTORY}" COMPONENT doc)
  ENDIF(fpath)
endforeach()















