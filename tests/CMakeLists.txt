# Test file list
set(TEST_SOURCES
    AESolutionSetTest.cpp
    DebugTest2.cpp
    EigenTest.cpp
    ExpandingMapTest.cpp
    FixedIntervalTest.cpp
    IntervalVectorTest.cpp
    KaucherArithmeticTest.cpp
    OriginalKrawczykTest.cpp
    SimpleTest.cpp
    ToleranceEmbeddingCorrect.cpp
    ToleranceEmbeddingTest.cpp
)

# Create a single test executable
add_executable(AllTests ${TEST_SOURCES})

# Link with Google Test and IntervalKrawczyk library
target_link_libraries(AllTests PRIVATE 
    IntervalKrawczyk 
    GTest::gtest_main
    GTest::gmock
)

# Enable CTest
enable_testing()

# Add test to CTest
add_test(NAME AllTests COMMAND AllTests)
