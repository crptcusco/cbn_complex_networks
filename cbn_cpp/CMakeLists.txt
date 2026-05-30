cmake_minimum_required(VERSION 3.10)
project(cbn_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include)

file(GLOB SOURCES "src/*.cpp")

add_library(cbn_core ${SOURCES})

add_executable(verify_cbn examples/verify_cbn.cpp)
target_link_libraries(verify_cbn cbn_core)

add_executable(run_pipeline examples/run_pipeline.cpp)
target_link_libraries(run_pipeline cbn_core)

add_executable(experiment_runner examples/experiment_runner.cpp)
target_link_libraries(experiment_runner cbn_core)
