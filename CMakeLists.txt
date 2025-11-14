cmake_minimum_required(VERSION 3.20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(
    cool-serial
    VERSION 0.0.0
    LANGUAGES CXX)
include_directories(${PROJECT_SOURCE_DIR}/include)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/v1.17.0.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

file (GLOB_RECURSE TEST_SRC ${PROJECT_SOURCE_DIR}/test/*.cpp)

add_executable(
    coolSerialTest
    ${TEST_SRC}
)

target_link_libraries(
    coolSerialTest
    GTest::gtest_main
    GTest::gmock
)

include (GoogleTest test)
gtest_discover_tests(coolSerialTest)
