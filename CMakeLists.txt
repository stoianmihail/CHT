cmake_minimum_required(VERSION 3.10)
project(cht)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -Wall -Wextra")

find_package(Threads REQUIRED)
set(THREADS_PREFER_PTHREAD_FLAG ON)

file(GLOB INCLUDE_H "include/cht/*.h")
set(EXAMPLE_FILES example.cc)

add_executable(example ${INCLUDE_H} ${EXAMPLE_FILES})
