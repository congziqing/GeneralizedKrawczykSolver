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

# Create executable for each test file
foreach(TEST_SOURCE ${TEST_SOURCES})
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    target_link_libraries(${TEST_NAME} PRIVATE IntervalKrawczyk)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
