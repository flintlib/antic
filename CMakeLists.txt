cmake_minimum_required(VERSION 3.1)

project(antic C)

file(READ "${CMAKE_CURRENT_SOURCE_DIR}/configure" CONFIGURE_CONTENTS)
string(REGEX MATCH "ANTIC_MAJOR=([0-9]*)" _ ${CONFIGURE_CONTENTS})
set(ANTIC_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "ANTIC_MINOR=([0-9]*)" _ ${CONFIGURE_CONTENTS})
set(ANTIC_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "ANTIC_PATCH=([0-9]*)" _ ${CONFIGURE_CONTENTS})
set(ANTIC_PATCH ${CMAKE_MATCH_1})

set(ANTIC_VERSION "${ANTIC_MAJOR}.${ANTIC_MINOR}.${ANTIC_PATCH}")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type (Debug, Release)" FORCE)
endif ()

set (BUILD_SHARED_LIBS yes CACHE BOOL "Build shared library or not")
set (BUILD_TESTING no CACHE BOOL "Build tests or not")

if (NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
    message("${CMAKE_BUILD_TYPE}")
    message(FATAL_ERROR "CMAKE_BUILD_TYPE must be one of: Debug, Release (current value: '${CMAKE_BUILD_TYPE}')")
endif ()

if ("${CMAKE_SIZEOF_VOID_P}" STREQUAL "8")
    set (PLATFORM "x64")
else ()
    set (PLATFORM "Win32")
endif()

if (BUILD_SHARED_LIBS)
    set (LIBRARY_TYPE dll)
else ()
    set (LIBRARY_TYPE lib)
endif()

set (MSVC_USE_MT yes CACHE BOOL "Use MT flags when compiling in MSVC")
set (MSVC_WARNING_LEVEL 1 CACHE STRING "MSVC warning level")

if (MSVC)
    if (MSVC_USE_MT)
        foreach(CompilerFlag CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE)
          string(REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
        endforeach()
    endif()
    foreach(CompilerFlag CMAKE_C_FLAGS)
      set(${CompilerFlag} "${${CompilerFlag}} /W${MSVC_WARNING_LEVEL}")
    endforeach()

    foreach(CompilerFlag CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE)
        if (BUILD_SHARED_LIBS)
            set(${CompilerFlag} "${${CompilerFlag}} /DMSC_USE_DLL")
        else ()
            set(${CompilerFlag} "${${CompilerFlag}} /DPTW32_STATIC_LIB")
        endif()
    endforeach()

    if ("${LIBRARY_TYPE}" STREQUAL "dll")
        set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS yes)
    endif()
endif()

set (DEPS flint2 mpfr mpir)

set (mpir_lib gmp mpir)
set (mpfr_lib mpfr)
set (flint2_lib ${LIBRARY_TYPE}_flint flint)

set (mpir_header gmp.h)
set (mpfr_header mpfr.h)
set (flint2_header flint/flint.h)

foreach (LIB ${DEPS})
    string (TOUPPER ${LIB} LIB_UPPER)
    find_library(${LIB_UPPER}_LIBRARY NAMES ${${LIB}_lib} HINTS ../${LIB}
                 PATH_SUFFIXES ${LIBRARY_TYPE}/${PLATFORM}/${CMAKE_BUILD_TYPE})
    if (NOT ${LIB_UPPER}_LIBRARY)
        message(FATAL_ERROR "${LIB} library not found.")
    endif()
    add_library(${LIB} UNKNOWN IMPORTED)
    set_property(TARGET ${LIB} PROPERTY IMPORTED_LOCATION ${${LIB_UPPER}_LIBRARY})
    message("${LIB} found in ${${LIB_UPPER}_LIBRARY}")
endforeach ()

foreach (LIB ${DEPS})
    string(TOUPPER ${LIB} HEADER_PKG)
    set (HEADER ${${LIB}_header})
    find_path(${HEADER_PKG}_INCLUDE_DIR NAMES ${HEADER} HINTS ../${LIB}
             PATH_SUFFIXES ${LIBRARY_TYPE}/${PLATFORM}/${CMAKE_BUILD_TYPE})
    if (NOT ${HEADER_PKG}_INCLUDE_DIR)
        message(FATAL_ERROR "${HEADER} header not found.")
    endif()
    message("${HEADER} found in ${${HEADER_PKG}_INCLUDE_DIR}")
    set (DEP_INCLUDE_DIRS ${DEP_INCLUDE_DIRS} ${${HEADER_PKG}_INCLUDE_DIR})
endforeach ()

file(GLOB TEMP "*.h")

foreach (TEMP_H ${TEMP})
    get_filename_component(FOLDER ${TEMP_H} NAME_WE)
    set(FOLDERS ${FOLDERS} ${FOLDER})
endforeach()

foreach (FOLDER ${FOLDERS})
    file(GLOB TEMP "${FOLDER}/*.c")
    set(SRC ${SRC} ${TEMP})
endforeach ()

include_directories(BEFORE ${antic_SOURCE_DIR})
include_directories(BEFORE ${DEP_INCLUDE_DIRS})

add_library(antic ${SRC})
target_link_libraries(antic ${DEPS})
target_compile_definitions(antic PRIVATE "ANTIC_BUILD_DLL")

set_target_properties(antic PROPERTIES VERSION ${ANTIC_VERSION} SOVERSION ${ANTIC_MAJOR})
if(WIN32)
    set_target_properties(antic PROPERTIES RUNTIME_OUTPUT_NAME "antic-${ANTIC_MAJOR}")
endif()

if(NOT MSVC)
    target_link_libraries(antic m)
endif()

install(TARGETS antic
            RUNTIME DESTINATION bin
            ARCHIVE DESTINATION lib
            LIBRARY DESTINATION lib
        )

foreach (FOLDER ${FOLDERS})
    set(HEADERS ${HEADERS} ${FOLDER}.h)
endforeach ()

install(FILES ${HEADERS} DESTINATION include/antic)

if (BUILD_TESTING)
    enable_testing()
    foreach (FOLDER ${FOLDERS})
        file(GLOB TEMP "${FOLDER}/test/*.c")
        foreach (TEST_SOURCE ${TEMP})
            get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
            add_executable(${FOLDER}-${TEST_NAME} ${TEST_SOURCE})
            target_link_libraries(${FOLDER}-${TEST_NAME} antic)
            add_test(${FOLDER}-${TEST_NAME} ${FOLDER}-${TEST_NAME})
        endforeach ()
    endforeach ()
endif ()
