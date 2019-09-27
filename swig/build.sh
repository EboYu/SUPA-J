#!/bin/bash

swig -c++ -I/../SVF/include/ -I/home/yinbo/LLVM/llvm-9.0.0.src/include -java SUPA.i
clang++ -fPIC SUPA_wrap.cxx -I$JAVA_HOME/include -I$JAVA_HOME/include/linux -I../SVF/include -I$LLVM_SRC/include -I$LLVM_DIR/include
