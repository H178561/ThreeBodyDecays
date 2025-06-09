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
    OrientationAngles.cpp
    ThreeBodyAmplitudeModel.cpp
    ClebschGordan.cpp
    #Examples/AmplitureModelExample.cpp
    #Examples/CombinedAmplitureModelExample.cpp
)

set_target_properties(ThreeBodyDecaysLib PROPERTIES POSITION_INDEPENDENT_CODE ON)

#target_link_libraries(ThreeBodyDecays ThreeBodyDecaysLib)

# Try to find GTest, but make it optional
find_package(GTest QUIET)
if(GTEST_FOUND)
    include_directories(${GTEST_INCLUDE_DIRS})

    # Add pthread for Google Test
    find_package(Threads REQUIRED)

    # Test executable
    add_executable(ThreeBodyDecaysTest
        Tests/ThreeBodyDecaysTest.cpp
        Tests/OrientationAnglesTest.cpp
        Tests/ThreeBodyAmplitudeModelTest.cpp
    )
    target_link_libraries(ThreeBodyDecaysTest
        ThreeBodyDecaysLib
        ${GTEST_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
    )

    # Add the CG test if GTest is found
    add_executable(ClebschGordanTest Tests/ClebschGordanTest.cpp)
    add_executable(JsonModelTest Tests/JsonModelTest.cpp)

    target_link_libraries(ClebschGordanTest ${GTEST_LIBRARIES} pthread ThreeBodyDecaysLib)
    add_test(NAME ClebschGordanTest COMMAND ClebschGordanTest)

    target_link_libraries(JsonModelTest ${GTEST_LIBRARIES} pthread ThreeBodyDecaysLib)
    add_test(NAME JsonModelTest COMMAND JsonModelTest)
    # Enable testing
    enable_testing()
    add_test(NAME ThreeBodyDecaysTest COMMAND ThreeBodyDecaysTest)

    message(STATUS "Google Test found, test target will be built")
else()
    message(STATUS "Google Test not found, test target will NOT be built")
endif()

# Füge das Beispiel hinzu
add_executable(AmplitudeModelExample Examples/AmplitudeModelExample.cpp)
target_link_libraries(AmplitudeModelExample ThreeBodyDecaysLib)
