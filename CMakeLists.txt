cmake_minimum_required(VERSION 3.10)
project(Epi)

# Compilation step for static resources
FUNCTION(ADD_RESOURCES out_var)
    SET(result)
    FOREACH(in_f ${ARGN})
        FILE(RELATIVE_PATH src_f ${CMAKE_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/${in_f})
        SET(out_f "${PROJECT_BINARY_DIR}/${in_f}.o")
        ADD_CUSTOM_COMMAND(OUTPUT ${out_f}
                COMMAND ld -r -b binary -o ${out_f} ${src_f}
                DEPENDS ${in_f}
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMENT "Building resource object ${out_f}"
                VERBATIM
                )
        LIST(APPEND result ${out_f})
    ENDFOREACH()
    SET(${out_var} "${result}" PARENT_SCOPE)
ENDFUNCTION()

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "-Wall -Wextra -O3")
SET(Boost_USE_STATIC_LIBS ON)

FIND_PACKAGE(Boost REQUIRED COMPONENTS program_options system log iostreams thread)
FIND_PACKAGE (Threads)

FIND_LIBRARY(ZLIB z)

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(include)

LINK_DIRECTORIES(${BOOST_LIBRARY_DIRS})

FILE(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp")
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/resources)

ADD_RESOURCES(RESOURCES resources/rfBm.txt.gz)
ADD_EXECUTABLE(Epi ${SOURCES} ${RESOURCES})

TARGET_LINK_LIBRARIES(Epi ${Boost_LIBRARIES} Threads::Threads z stdc++fs)