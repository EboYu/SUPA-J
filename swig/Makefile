
.PHONY: clean wrap compile build 

clean:
	rm -rf *.java *.cxx

wrap:
	rm -rf *.java *.cxx 
	swig -c++ -I/../SVF/include/ -I${LLVM_SRC}/include -java SUPA.i

compile:
	clang++ -fPIC SUPA_wrap.cxx -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux -I/home/yinbo/disk/workspace/SUPA-J/SVF/include -I${LLVM_SRC}/include -I${LLVM_DIR}/include

build:
	wrap
	compile


