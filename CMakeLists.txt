
set(the_description "The Core Functionality")

ut_add_module(utcore)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(tracing_hdr_files "")
set(tracing_src_files "")

set(tracing_extra_include_dirs "")
set(tracing_extra_libraries "")


# @todo are tracing headers installed correctly ??
IF(ENABLE_TRACING_DTRACE)
    set(tracing_hdr_files ${tracing_hdr_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_dtrace.h")
ENDIF(ENABLE_TRACING_DTRACE)

IF(ENABLE_TRACING_ETW)
    set(tracing_hdr_files ${tracing_hdr_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.h")
    set(tracing_src_files ${tracing_src_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.rc")
ENDIF(ENABLE_TRACING_ETW)

IF(ENABLE_TRACING_LTTNGUST)
    set(tracing_extra_include_dirs ${tracing_extra_include_dirs} ${LTTNGUST_INCLUDE_DIRS})
    set(tracing_extra_libraries LTTng::UST)
ENDIF(ENABLE_TRACING_LTTNGUST)


ut_module_include_directories(${UBITRACK_CORE_DEPS_INCLUDE_DIR} ${tracing_extra_include_dirs})
ut_glob_module_sources(HEADERS "src/*.h" "src/*/*.h" "src/*/*/*.h" "src/*/*/*/*.h" "src/*/*/*/*/*.h" ${tracing_hdr_files} SOURCES "src/*/*.cpp" "src/*/*/*.cpp" "src/*/*/*/*.cpp" "src/*/*/*/*/*.cpp" ${tracing_src_files})
ut_create_module(${TINYXML_LIBRARIES} ${LOG4CPP_LIBRARIES} ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} ${MSGPACK_LIBRARIES} ${tracing_extra_libraries})

ut_add_module_tests()

IF(ENABLE_TRACING_DTRACE)
    DTRACE_INSTRUMENT(utcore)
ENDIF(ENABLE_TRACING_DTRACE)

IF(ENABLE_TRACING_ETW)
    ETW_INSTRUMENT(utcore)
ENDIF(ENABLE_TRACING_ETW)


# install custom files
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