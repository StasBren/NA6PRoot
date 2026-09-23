cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(mwpc_phase_a LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT TARGET Garfield::Garfield)
  find_package(Garfield REQUIRED)
endif()

add_executable(mwpc_phase_a mwpc_phase_a.cxx)
target_link_libraries(mwpc_phase_a PRIVATE Garfield::Garfield)
