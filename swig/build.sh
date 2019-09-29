#!/bin/bash

swig -c++ -I/../SVF/include/ -I${LLVM_SRC}/include -java SUPA.i
clang++ -fPIC SUPA_wrap.cxx -I$JAVA_HOME/include -I$JAVA_HOME/include/linux -I../SVF/include -I$LLVM_SRC/include -I$LLVM_DIR/include
