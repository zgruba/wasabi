cmake_minimum_required(VERSION 3.12)
project(WasabiProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build (Debug or Release)" FORCE)
endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)   # static libs
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)   # shared libs
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)   # executables

set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -fsanitize=address,undefined")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

# -------------------------------
# Dependencies
# -------------------------------

# Replace with your real paths:
set(ISOSPEC_INCLUDE_DIR "/home/zosia/Studies/wasabi_dir/IsoSpec/")
set(ISOSPEC_LIB_DIR "/home/zosia/Studies/wasabi_dir/IsoSpec/IsoSpec++/")

include_directories(${ISOSPEC_INCLUDE_DIR})
link_directories(${ISOSPEC_LIB_DIR})

find_package(OpenMP REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options log)
# OpenMs has some libsvm conflicts :)
#find_package(OpenMS REQUIRED)

# -------------------------------
# Libraries
# -------------------------------
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${Boost_INCLUDE_DIR})

# --- Header-only amino acid table ---
add_library(amino_acids INTERFACE)
target_include_directories(amino_acids INTERFACE ${CMAKE_SOURCE_DIR}/include)

# --- Header-only PTM library ---
add_library(ptm INTERFACE)
target_include_directories(ptm INTERFACE ${CMAKE_SOURCE_DIR}/include)

# --- Wasabi DB Templates  ---
add_library(wasabi_db INTERFACE)
target_include_directories(wasabi_db INTERFACE ${CMAKE_SOURCE_DIR}/include)

# --- IsoSpec ---
add_library(iso_lib
        source/iso.cpp
        include/iso.hpp
)
target_link_libraries(iso_lib PRIVATE IsoSpec++ OpenMP::OpenMP_CXX amino_acids)


# --- OpenMS ---
add_library(openms_lib
        source/parsing.cpp
        include/parsing.hpp
)
target_link_libraries(openms_lib PRIVATE OpenMS OpenMP::OpenMP_CXX amino_acids)

# -------------------------------
# Executables
# -------------------------------

add_executable(prepare_csv_for_deeplc source/prepare_csv_for_deeplc.cpp)
target_link_libraries(prepare_csv_for_deeplc PRIVATE iso_lib openms_lib Boost::program_options)

add_executable(wasabi_annotate source/wasabi_annotate.cpp)
target_link_libraries(wasabi_annotate PRIVATE iso_lib IsoSpec++ openms_lib OpenMS wasabi_db OpenMP::OpenMP_CXX amino_acids Boost::program_options Boost::log)
