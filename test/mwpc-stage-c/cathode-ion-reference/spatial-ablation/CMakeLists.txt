# Standalone build of the S-ladder C++ tools against the FROZEN sources on the wire branch (unchanged).
#   cmake -S . -B build -DFROZEN_DIR=<wire-branch>/test/mwpc-stage-c/cathode-ion-reference && cmake --build build
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(NA6PSpatialAblation LANGUAGES CXX)
find_package(Garfield REQUIRED)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
set(FROZEN_DIR "" CACHE PATH "directory with p3_export_full_spectral_paths.cc, Plane2100Archive.hh, Plane2100Response.hh")
if(NOT EXISTS "${FROZEN_DIR}/p3_export_full_spectral_paths.cc")
  message(FATAL_ERROR "set -DFROZEN_DIR to the wire branch's test/mwpc-stage-c/cathode-ion-reference")
endif()
file(SHA256 "${FROZEN_DIR}/p3_export_full_spectral_paths.cc" FROZEN_SHA)
message(STATUS "frozen exporter source sha256 ${FROZEN_SHA}")
add_executable(p3_green_grid_paths p3_green_grid_paths.cc)
target_include_directories(p3_green_grid_paths PRIVATE "${FROZEN_DIR}")
target_link_libraries(p3_green_grid_paths PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(p3_green_grid_paths PRIVATE cxx_std_17)
target_compile_definitions(p3_green_grid_paths PRIVATE FROZEN_EXPORTER_SHA256="${FROZEN_SHA}")

# Birth-record -> K2A cloud moments (no ion transport). Reads archives with the frozen Plane2100Archive.hh reader.
add_executable(p3_birth_moments p3_birth_moments.cc)
target_include_directories(p3_birth_moments PRIVATE "${FROZEN_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(p3_birth_moments PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
target_compile_features(p3_birth_moments PRIVATE cxx_std_17)
