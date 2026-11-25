# K2A inline full-birth moments (milestone 1): separate CMake project; no existing CMake file is changed.
#   cmake -S <this dir> -B build [-DCMAKE_BUILD_TYPE=<same as your frozen generator build>] && cmake --build build -j2
#   -DK2A_INLINE_T0_ONLY=ON builds only the dependency-free T0 unit test (no Garfield++, no ROOT).
# Targets:
#   mwpc_p3_noFW_penning057_k2a_inline_v1  copy of the frozen 2100 V generator + marked K2A-INLINE blocks
#                                          (linked exactly like the frozen target: Garfield::Garfield ROOT::Core/RIO/Tree)
#   k2a_archive_dump                       read-only TTree dump of an archive (independent T1 reference input)
#   p3_birth_moments                       spatial-ablation batch moments (frozen reader + BirthMoments.hh), full or weighted
#   moment_accumulator_t0                  T0 synthetic unit test
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(NA6PK2AInline LANGUAGES CXX)
option(K2A_INLINE_T0_ONLY "build only the T0 unit test (no Garfield++/ROOT)" OFF)
set(AV "${CMAKE_CURRENT_SOURCE_DIR}/../../../mwpc-avalanche")
set(CIR "${CMAKE_CURRENT_SOURCE_DIR}/..")
set(FROZEN_GEN "${AV}/phaseB_native_muon_p3_noFW_penning057_cathode_ion_reference_v1.cxx")
set(INLINE_GEN "${AV}/phaseB_native_muon_p3_noFW_penning057_k2a_inline_v1.cxx")
file(SHA256 "${FROZEN_GEN}" FROZEN_SHA)
if(NOT FROZEN_SHA STREQUAL "1436a71575cc34c002b49130a017409364e551ed11a35d1ee31df17224f5874f")
  message(FATAL_ERROR "frozen generator source sha256 ${FROZEN_SHA} differs from the recorded frozen source")
endif()

add_executable(moment_accumulator_t0 tests/moment_accumulator_t0.cc)
target_compile_features(moment_accumulator_t0 PRIVATE cxx_std_17)

if(NOT K2A_INLINE_T0_ONLY)
  find_package(Garfield REQUIRED)
  find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
  add_executable(mwpc_p3_noFW_penning057_k2a_inline_v1 "${INLINE_GEN}")
  target_link_libraries(mwpc_p3_noFW_penning057_k2a_inline_v1 PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
  target_compile_features(mwpc_p3_noFW_penning057_k2a_inline_v1 PRIVATE cxx_std_17)

  add_executable(k2a_archive_dump k2a_archive_dump.cc)
  target_link_libraries(k2a_archive_dump PRIVATE ROOT::Core ROOT::RIO ROOT::Tree)
  target_compile_features(k2a_archive_dump PRIVATE cxx_std_17)

  add_executable(p3_birth_moments "${CIR}/spatial-ablation/p3_birth_moments.cc")
  target_include_directories(p3_birth_moments PRIVATE "${CIR}" "${CIR}/spatial-ablation")
  target_link_libraries(p3_birth_moments PRIVATE Garfield::Garfield ROOT::Core ROOT::RIO ROOT::Tree)
  target_compile_features(p3_birth_moments PRIVATE cxx_std_17)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" BT)
foreach(f IN ITEMS "${INLINE_GEN}" "${CMAKE_CURRENT_SOURCE_DIR}/MomentAccumulator.hh" "${CIR}/CathodeIonArchive.hh"
                   "${CIR}/Plane2100Archive.hh" "${CIR}/spatial-ablation/BirthMoments.hh" "${CIR}/spatial-ablation/p3_birth_moments.cc")
  file(SHA256 "${f}" s)
  get_filename_component(n "${f}" NAME)
  string(APPEND SRC_SHAS "source_sha256 ${n} ${s}\n")
endforeach()
file(WRITE "${CMAKE_BINARY_DIR}/k2a_inline_build_info.txt"
  "compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} ${CMAKE_CXX_COMPILER}\n"
  "build_type '${CMAKE_BUILD_TYPE}'\ncxx_flags '${CMAKE_CXX_FLAGS}'\ncxx_flags_build_type '${CMAKE_CXX_FLAGS_${BT}}'\n"
  "t0_only ${K2A_INLINE_T0_ONLY}\nfrozen_generator_sha256 ${FROZEN_SHA}\n${SRC_SHAS}")
message(STATUS "K2A inline build info: ${CMAKE_BINARY_DIR}/k2a_inline_build_info.txt (build type '${CMAKE_BUILD_TYPE}')")
