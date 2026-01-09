# 演示程序列表
set(DEMO_SOURCES
    ComplexPolynomialDemo.cpp
    FixedDemo.cpp
    GeneralizedKrawczykAdvancedDemo.cpp
    GeneralizedKrawczykDemo.cpp
    GeneralizedKrawczykHomotopy.cpp
    HomotopyDemoSimple.cpp
    KrawczykDemo.cpp
    ToleranceEmbeddingFinal.cpp
    ToleranceEmbeddingFixed.cpp
    ToleranceEmbeddingFixedSimple.cpp
    ToleranceEmbeddingHomotopy.cpp
)

# 为每个演示文件创建可执行文件
foreach(DEMO_SOURCE ${DEMO_SOURCES})
    get_filename_component(DEMO_NAME ${DEMO_SOURCE} NAME_WE)
    add_executable(${DEMO_NAME} ${DEMO_SOURCE})
    target_link_libraries(${DEMO_NAME} PRIVATE IntervalKrawczyk)
endforeach()
