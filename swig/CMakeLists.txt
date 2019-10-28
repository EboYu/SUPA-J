cmake_minimum_required(VERSION 3.4.3)

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE($(SWIG_USE_FILE))

INCLUDE_DIRECTORIES(${JAVA_HOME}/include)

INCLUDE_DIRECTORIES(${JAVA_HOME}/include/linux)
set_property(SOURCE 
    include/DDA/DDAClient.i
    include/DDA/DDAPass.i
    include/DDA/DDAStat.i
    include/DDA/DDAVFSolver.i
    include/MemoryModel/PointerAnalysis.i
    include/Util/BasicTypes.i
    include/Util/PTAStat.i
    include/Util/SVFModule.i
    include/Util/SVFUtil.i
    PROPERTY CPLUSPLUS ON
    )
swig_add_library(libSUPA LANGUAGE java SOURCEs 
    include/DDA/DDAClient.i
    include/DDA/DDAPass.i
    include/DDA/DDAStat.i
    include/DDA/DDAVFSolver.i
    include/MemoryModel/PointerAnalysis.i
    include/Util/BasicTypes.i
    include/Util/PTAStat.i
    include/Util/SVFModule.i
    include/Util/SVFUtil.i)


link_directories( ${CMAKE_BINARY_DIR}/lib/SUPA )
target_link_libraries(libSUPA)