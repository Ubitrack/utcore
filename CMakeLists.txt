set(the_description "The Core Functionality")

ut_add_module(utcore)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
UT_OPTION(UT_ENABLE_DTRACE "Enable DTrace tracing of eventqueue activity?" OFF)
UT_OPTION(UT_ENABLE_ETW "Enable ETW Tracing of eventqueue activity?" OFF)
UT_OPTION(UT_ENABLE_LTTNGUST "Enable LTTNG-UST Tracing of eventqueue activity?" OFF)

set(tracing_hdr_files "")
set(tracing_src_files "")

set(tracing_extra_include_dirs "")
set(tracing_extra_libraries "")


IF(UT_ENABLE_DTRACE)
    MESSAGE( STATUS "Enabled DTrace.")
    INCLUDE(dtrace)
    add_definitions("-DHAVE_DTRACE")
    add_definitions("-DENABLE_EVENT_TRACING")
    include_directories(${CMAKE_BINARY_DIR})
    set(tracing_hdr_files ${tracing_hdr_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_dtrace.h")
ENDIF(UT_ENABLE_DTRACE)

IF(UT_ENABLE_LTTNGUST)
    MESSAGE( STATUS "Enabled LTTNG-UST.")
    find_package(LTTngUST)
    if(LTTNGUST_FOUND)
        add_definitions("-DHAVE_LTTNGUST")
        add_definitions("-DENABLE_EVENT_TRACING")
        include_directories(${CMAKE_BINARY_DIR})
        set(tracing_extra_include_dirs ${tracing_extra_include_dirs} ${LTTNGUST_INCLUDE_DIRS})
        set(tracing_extra_libraries ${tracing_extra_libraries} ${LTTNGUST_LIBRARIES})
    endif(LTTNGUST_FOUND)
ENDIF(UT_ENABLE_LTTNGUST)


IF(UT_ENABLE_ETW)
    MESSAGE( STATUS "Enabled Event Tracing for Windows.")
    INCLUDE(etw)
    add_definitions("-DHAVE_ETW")
    add_definitions("-DENABLE_EVENT_TRACING")
    include_directories(${CMAKE_BINARY_DIR})
    set(tracing_hdr_files ${tracing_hdr_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.h")
    set(tracing_src_files ${tracing_src_files} "${CMAKE_BINARY_DIR}/utUtil/probes_ubitrack_etw.rc")
ENDIF(UT_ENABLE_ETW)


ut_module_include_directories(${UBITRACK_CORE_DEPS_INCLUDE_DIR} ${tracing_extra_include_dirs})
ut_glob_module_sources(HEADERS "src/*.h" "src/*/*.h" "src/*/*/*.h" "src/*/*/*/*.h" "src/*/*/*/*/*.h" SOURCES "src/*/*.cpp" "src/*/*/*.cpp" "src/*/*/*/*.cpp" "src/*/*/*/*/*.cpp")
ut_create_module(${TINYXML_LIBRARIES} ${LOG4CPP_LIBRARIES} ${LAPACK_LIBRARIES} ${Boost_LIBRARIES} ${tracing_extra_libraries})

ut_add_module_tests()

IF(UT_ENABLE_DTRACE)
    DTRACE_INSTRUMENT(utcore)
ENDIF(UT_ENABLE_DTRACE)

# do we need to instrument lttng probes ??

IF(UT_ENABLE_ETW)
    ETW_INSTRUMENT(utcore)
ENDIF(UT_ENABLE_ETW)
