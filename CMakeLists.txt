cmake_minimum_required(VERSION 3.10)
project(ThreeBodyDecays)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Add source files
add_library(ThreeBodyDecaysLib
    src/ThreeBodyDecays.cpp
    src/OrientationAngles.cpp
    src/ThreeBodyAmplitudeModel.cpp
    src/ClebschGordan.cpp
    src/FormFactors.cpp
    src/WignerFunctions.cpp
    src/SpinParity.cpp
    src/WignerRotation.cpp
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

    # Define test source files
    set(TEST_SOURCES
        Tests/ThreeBodyDecaysTest.cpp
        Tests/OrientationAnglesTest.cpp
        Tests/ThreeBodyAmplitudeModelTest.cpp
        Tests/ClebschGordanTest.cpp
        Tests/JsonModelTest.cpp
        Tests/BWTest.cpp
        Tests/FormFactorsTest.cpp
        Tests/WignerFunctionsTest.cpp
    )

    # Individual test executables
    add_executable(ThreeBodyDecaysTest Tests/ThreeBodyDecaysTest.cpp)
    target_link_libraries(ThreeBodyDecaysTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(ClebschGordanTest Tests/ClebschGordanTest.cpp)
    target_link_libraries(ClebschGordanTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(JsonModelTest Tests/JsonModelTest.cpp)
    target_link_libraries(JsonModelTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(BWTest Tests/BWTest.cpp)
    target_link_libraries(BWTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(FormFactorsTest Tests/FormFactorsTest.cpp)
    target_link_libraries(FormFactorsTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(OrientationAnglesTest Tests/OrientationAnglesTest.cpp)
    target_link_libraries(OrientationAnglesTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(ThreeBodyAmplitudeModelTest Tests/ThreeBodyAmplitudeModelTest.cpp)
    target_link_libraries(ThreeBodyAmplitudeModelTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    add_executable(WignerFunctionsTest Tests/WignerFunctionsTest.cpp)
    target_link_libraries(WignerFunctionsTest ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    # Combined test executable with a dedicated main function
    add_executable(RunTests Tests/test_main.cpp ${TEST_SOURCES})
    target_compile_definitions(RunTests PRIVATE RUNNING_COMBINED_TESTS)
    target_link_libraries(RunTests ThreeBodyDecaysLib ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

    # Enable testing
    enable_testing()
    add_test(NAME ThreeBodyDecaysTest COMMAND ThreeBodyDecaysTest)
    add_test(NAME ClebschGordanTest COMMAND ClebschGordanTest)
    add_test(NAME JsonModelTest COMMAND JsonModelTest)
    add_test(NAME BWTest COMMAND BWTest)
    add_test(NAME FormFactorsTest COMMAND FormFactorsTest)
    add_test(NAME BlattWeisskopfTest COMMAND BlattWeisskopfTest)
    add_test(NAME OrientationAnglesTest COMMAND OrientationAnglesTest)
    add_test(NAME ThreeBodyAmplitudeModelTest COMMAND ThreeBodyAmplitudeModelTest)
    add_test(NAME WignerFunctionsTest COMMAND WignerFunctionsTest)
    add_test(NAME RunAllTests COMMAND RunTests)

    message(STATUS "Google Test found, test targets will be built")
else()
    message(STATUS "Google Test not found, test targets will NOT be built")
endif()

# Füge das Beispiel hinzu
add_executable(AmplitudeModelExample Examples/AmplitudeModelExample.cpp)
target_link_libraries(AmplitudeModelExample ThreeBodyDecaysLib)
