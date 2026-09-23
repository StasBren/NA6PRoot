cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(NA6PMWPCAvalanchePhaseA LANGUAGES CXX)

find_package(Garfield REQUIRED)

add_executable(mwpc_phase_a phaseA_single_electron.cxx)
target_link_libraries(mwpc_phase_a PRIVATE Garfield::Garfield)
target_compile_features(mwpc_phase_a PRIVATE cxx_std_17)
