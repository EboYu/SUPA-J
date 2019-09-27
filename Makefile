# Makefile for SUPA-J

# Path to SUPA such that headers for the C bindings exist
SUPA_HOME=SVF/include
JNAERATOR_VERSION=0.12-shaded
JAVA=java

.PHONY: help clean binding

.DEFAULT_GOAL := help
	
help:
	@echo "Use \`make <target> [SUPA_HOME=<path>]' where <path> is the"
	@echo "location of SUPA such that headers for the C-bindings"
	@echo "exist in \$$SUPA_HOME/include and <target> is one of"
	@echo "    help     to display this help message"
	@echo "    clean    to remove existing bindings"
	@echo "    binding     to generate Java bindings to SUPA using JNAerator"

clean:
	-rm src/main/java/org/supa/binding/SUPALibrary.java

binding: clean src/main/java/org/supa/binding/SUPALibrary.java

src/main/java/org/supa/binding/SUPALibrary.java: jnaerator.jar config.jnaerator
	$(JAVA) -jar $<
	sed -i 's/@Library("SUPA")/@Library("SUPA")/' $@

jnaerator.jar:
ifeq ("$(wildcard jnaerator.jar)","")
	wget http://jnaerator.googlecode.com/files/jnaerator-$(JNAERATOR_VERSION).jar -O $@
endif
