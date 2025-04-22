cmake_minimum_required(VERSION 3.10)
project(ThreeBodyDecays)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add source files
add_library(ThreeBodyDecaysLib
    ThreeBodyDecays.cpp
)

target_link_libraries(ThreeBodyDecays ThreeBodyDecaysLib)

# Try to find GTest, but make it optional
find_package(GTest QUIET)
if(GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})

    # Add pthread for Google Test
    find_package(Threads REQUIRED)

    # Test executable
    add_executable(ThreeBodyDecaysTest
        Tests/ThreeBodyDecaysTest.cpp
    )
    target_link_libraries(ThreeBodyDecaysTest
        ThreeBodyDecaysLib
        ${GTEST_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
    )

    # Enable testing
    enable_testing()
    add_test(NAME ThreeBodyDecaysTest COMMAND ThreeBodyDecaysTest)

    message(STATUS "Google Test found, test target will be built")
else()
    message(STATUS "Google Test not found, test target will NOT be built")
endif()
