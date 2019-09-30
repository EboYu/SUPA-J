#!/usr/bin/python
import os
import subprocess
import sys
import shutil

standard_include = '/usr/include/'

pointerType = {'PointerAnalysis','Function','GlobalVariable', 'GlobalAlias', 'LLVMContext', 'Module', 'LLVMModuleSet', 'Instruction', 'Value', 'CallSite',
			   'PointerType', 'BasicBlock', 'ConstantExpr', 'u32_t','PTACallGraph','PAG','PTAStat','CallGraphSCC','ICFG','CHGraph','TypeSystem','MemObj',
			   'CallInst','PTACallGraphNode','PTDataTy','DiffPTDataTy','IncDFPTDataTy','SVFG','SVFGSCC','PAGNode','SVFGEdge','AndersenWaveDiff','DDAStat',
			   'IndirectSVFGEdge','DirectSVFGEdge','LoadSVFGNode','StoreSVFGNode','AllocaInst','GepSVFGNode','AddrSVFGNode','DDAClient','DominatorTree'}

def header_include(name,path,std):
	if std:
		return '<' + os.path.join(path[len(standard_include):],name + '.h') + '>'
	else:
		return '"' + path + '/' + name + '.h"'
		
def runprog(prog,args):
	popen_args = prog
	for arg in args:
		popen_args = popen_args + ' ' + arg
	proc = subprocess.Popen(popen_args,stdout=subprocess.PIPE,shell=True)
	while proc.poll() == None:
		pass
	return proc.stdout.read()
	
def join_paths(paths):
	def do_joins(path,paths):
		if len(paths) == 1:
			return os.path.join(path,paths[0])
		else:
			return do_joins(os.path.join(path,paths[0]),paths[1:])
	return do_joins(paths[0],paths[1:])
	
def last_token(string,delimiter):
	s = string.split(delimiter)
	if s[len(s)-1] != '':
		return s[len(s)-1]
	else:
		return s[len(s)-2]

def generate_swig_interfaces(header_path):
	assert(os.path.isdir(header_path))
	if header_path[len(header_path)-1] != '/':
		header_path = header_path + '/'
	std_include = False
	if header_path.find(standard_include) == 0:
		std_include = True
	headers = []
	for root, dirs, files in os.walk(header_path):
		for file in files:
			if root.find('DDA') >= 0 or root.find('WPA') >= 0 or file == 'BasicTypes.h' \
					or file == 'PointerAnalysis.h' or file == 'PTAStat.h' or file == 'SVFUtil.h'\
					or file == 'SVFModule.h':
				if file.rpartition('.')[2] == 'h':
					headers.append(os.path.join(root, file))
	print 'Generating JNI interface for headers in ' + header_path + ' using SWIG.'

	for header in headers:
		name = os.path.basename(header).rpartition('.')[0]
		relpath = os.path.normpath('./' + last_token(header_path,'/') + '/' + os.path.dirname(header)[len(header_path):] + '/')
		if relpath.find('/MTA')!=-1 or relpath.find('/sys')!=-1 or relpath.find('/bits')!=-1:
			continue
		print 'Preprocessing ' + header + ' for SWIG at ' + relpath
		if os.path.isdir(relpath) == False:
				os.makedirs(relpath)
		hfile = open(header)
		contents = hfile.read()
		hfile.close()
		out = open(join_paths([relpath,name + '.i']),"w")
		out.write('%module ' + name + '_SWIG\n%{\n#include ' + header_include(name,os.path.dirname(header),std_include) + '\n%}\n')

		if name == 'FlowDDA' or name == 'FlowSensitive' or name == 'Andersen':
			out.write('%import "../MemoryModel/PointerAnalysis.i"\n')
		if name == 'DDAStat' or name == 'WPAStat':
			out.write('%import "../Util/PTAStat.i"\n')
		if name == 'TypeAnalysis' or name == 'AndersenSFR':
			out.write('%import "Andersen.i"\n')

		if contents.find('public ModulePass') >=0:
			out.write('\n%import "llvm/Pass.h"\n')
			out.write('%import "../Util/BasicTypes.i"\n')

		out.write('%include <cpointer.i>\n')
		for type in pointerType:
			if contents.find(type+'*') >=0 or contents.find(type+' *') >=0:
				out.write('%pointer_class('+type+','+type+'P)\n')
		# if contents.find('Value*') >= 0:
		# 	out.write('%pointer_class(Value,ValueP)\n')
		# if contents.find('PointerAnalysis*') >= 0:
		# 	out.write('%pointer_class(PointerAnalysis,PointerAnalysisP)\n')
		# if contents.find('SVFG*') >= 0:
		# 	out.write('%pointer_class(SVFG,SVFGP)\n')
		# if contents.find('SVFGSCC*') >= 0:
		# 	out.write('%pointer_class(SVFGSCC,LSVFGSCCP)\n')
		# if contents.find('SVFGEdge*') >= 0:
		# 	out.write('%pointer_class(SVFGEdge,SVFGEdgeP)\n')
		# if contents.find('DDAClient*') >= 0:
		# 	out.write('%pointer_class(DDAClient,DDAClientP)\n')
		# if contents.find('DDAPass *') >=0:
		# 	out.write('%pointer_class(DDAPass,DDAPassP)\n')
		if contents.find('char **') >= 0:
			out.write('%include "carrays.i"\n')
			out.write('%array_functions(char *,StringArray)\n')
		# if contents.find('unsigned *') >= 0:
		# 	out.write('%pointer_class(unsigned,UnsignedIntArray)\n')
		out.write(contents)
		out.close()
		cwd = os.getcwd()
		packageName = relpath.replace("include/","")
	 	print 'swig ' + '-c++ -java -package org.supa.bindings.'+packageName+' -outdir src/org/supa/bindings/'+packageName+' -o jni/'+packageName+'/'+name+'_wrap.cxx ' + name + '.i'
		if os.path.isdir("src/org/supa/bindings/"+packageName) == False:
				os.makedirs("src/org/supa/bindings/"+packageName)
		if os.path.isdir("jni/"+packageName) == False:
				os.makedirs("jni/"+packageName)
		os.chdir(relpath)
		out = runprog("swig",["-I${LLVM_SRC}/include -I/home/yinbo/disk/workspace/SUPA-J/SVF/include/ -c++ -java -package org.supa.bindings."+packageName, "-outdir ../../src/org/supa/bindings/"+ packageName, "-o ../../jni/"+packageName+'/'+name+"_wrap.cxx",name + '.i', ])
		os.chdir(cwd)
		print out

def regenerate_swig_interfaces(header_path):
	if header_path[len(header_path)-1] != '/':
		header_path = header_path + '/'
	headers = []
	for root, dirs, files in os.walk(header_path):
		for file in files:
			if file.rpartition('.')[2] == 'i':
				headers.append(os.path.join(root, file))
	for header in headers:
		name = os.path.basename(header).rpartition('.')[0]
		relpath = os.path.normpath(
			'./' + last_token(header_path, '/') + '/' + os.path.dirname(header)[len(header_path):] + '/')
		cwd = os.getcwd()
		packageName = relpath.replace("include/", "")
		print 'swig ' + '-c++ -java -package org.supa.bindings.' + packageName + ' -outdir src/org/supa/bindings/' + packageName + ' -o jni/' + packageName + '/' + name + '_wrap.cxx ' + name + '.i'
		if os.path.isdir("src/org/supa/bindings/" + packageName) == False:
			os.makedirs("src/org/supa/bindings/" + packageName)
		if os.path.isdir("jni/" + packageName) == False:
			os.makedirs("jni/" + packageName)
		os.chdir(relpath)
		if name == 'SVFModule':
			out = runprog("swig", [
				"-I/home/yinbo/disk/workspace/SUPA-J/SVF/include/ -c++ -java -package org.supa.bindings." + packageName,
				"-outdir ../../src/org/supa/bindings/" + packageName,
				"-o ../../jni/" + packageName + '/' + name + "_wrap.cxx", name + '.i', ])
		else:
			out = runprog("swig", [
				"-I${LLVM_SRC}/include -I/home/yinbo/disk/workspace/SUPA-J/SVF/include/ -c++ -java -package org.supa.bindings." + packageName,
				"-outdir ../../src/org/supa/bindings/" + packageName,
				"-o ../../jni/" + packageName + '/' + name + "_wrap.cxx", name + '.i', ])
		os.chdir(cwd)
		print out

if __name__ == '__main__':
	# if os.path.isdir("include") == True:
	# 	shutil.rmtree("include/")
	# if os.path.isdir("src") == True:
	# 	shutil.rmtree("src/")
	# if os.path.isdir("jni") == True:
	# 	shutil.rmtree("jni/")
	# generate_swig_interfaces(sys.argv[1])
	# print 'Finish interface files generation'
	regenerate_swig_interfaces(sys.argv[2])
	print 'Generate wrap file'
