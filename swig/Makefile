
.PHONY: clean wrap compile build 
JAVA = "/opt/jdk1.8.0_191/bin/java"
JAVAC = "/opt/jdk1.8.0_191/bin/javac" -d .
JAVASRCS = "src/org/supa/bindings/*.java"

wrap:
	rm -rf src/ jni/
ifeq ("$(wildcard jni)", "")
		 mkdir jni
endif
ifeq ("$(wildcard src/org/supa/bindings)", "")
		 mkdir -p src/org/supa/bindings/
endif
	swig -c++ -I/home/yinbo/disk/workspace/SUPA-J/SVF/include/ -I${LLVM_SRC}/include -I/usr/include -I/usr/locale/include -java -package org.supa.bindings -outdir src/org/supa/bindings/ -fcompact -o jni/SUPA_wrap.cpp SUPA.i 

compile:
	clang++ -fPIC -c -ferror-limit=4 jni/SUPA_wrap.cpp -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux -I/home/yinbo/disk/workspace/SUPA-J/SVF/include -I${LLVM_DIR}/include -I${LLVM_SRC}/include 

java_compile:
ifeq ("$(wildcard target)", "")
		 mkdir target
endif
	$(JAVAC) $(JAVASRCS)

java_run:
	env LD_LIBRARY_PATH=$$PWD $(JAVA) runme

clean:
	-rm -rf *.java *.cxx

