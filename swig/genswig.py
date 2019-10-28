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

defined_headers = ['/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/BasicTypes.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/PTAStat.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/SVFUtil.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/SVFModule.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/MemoryModel/PointerAnalysis.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAVFSolver.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAStat.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAClient.h','/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAPass.h']

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
				"-I/home/yinbo/disk/workspace/SUPA-J/SVF/include/ -c++ -java -package org.supa.bindings.",
				"-outdir ../../src/org/supa/bindings/",
				"-o ../../jni/" + name + "_wrap.cxx", name + '.i', ])
		else:
			out = runprog("swig", [
				"-I${LLVM_SRC}/include -I/home/yinbo/disk/workspace/SUPA-J/SVF/include/ -c++ -java -package org.supa.bindings.",
				"-outdir ../../src/org/supa/bindings/",
				"-o ../../jni/" + name + "_wrap.cxx", name + '.i', ])
		os.chdir(cwd)
		print out

def generate_swig_single(header_path):
	assert (os.path.isdir(header_path))
	if header_path[len(header_path) - 1] != '/':
		header_path = header_path + '/'
	std_include = False
	if header_path.find(standard_include) == 0:
		std_include = True
	headers = []
	# for root, dirs, files in os.walk(header_path):
	# 	for file in files:
	# 		if root.find('DDA') >= 0 or file == 'BasicTypes.h' \
	# 				or file == 'PointerAnalysis.h' or file == 'PTAStat.h' or file == 'SVFUtil.h' \
	# 				or file == 'SVFModule.h':
	# 			if file.rpartition('.')[2] == 'h':
	# 				headers.append(os.path.join(root, file))
	print 'Generating JNI interface for headers in ' + header_path + ' using SWIG.'
	out = open('SUPA.i', "w")
	out.write('%module SUPA\n%{\n')
	content = str()
	sys_includes= set()
	includes = set()
	typedef = str()

	for header in defined_headers:
		# name = os.path.basename(header).rpartition('.')[0]
	# for header in headers:
		name = os.path.basename(header).rpartition('.')[0]
		relpath = os.path.normpath(
			'./' + last_token(header_path, '/') + '/' + os.path.dirname(header)[len(header_path):] + '/')
		if relpath.find('/MTA') != -1 or relpath.find('/sys') != -1 or relpath.find('/bits') != -1:
			continue
		print 'Preprocessing ' + header + ' for SWIG'
		hfile = open(header)

		while 1:
			line = hfile.readline()
			if not line:
				break
			if (line.startswith('#ifndef') and line.endswith('_H_\n')) or (line.startswith('#define') and line.endswith('_H_\n')) or line.startswith('#endif'):
				i=0
			elif line.startswith('#include'):
				if(line.startswith('#include <')):
					sys_includes.add(line)
				else:
					inName = line.replace('#include \"','').replace('\"\n','')
					length = len(headers)
					for h in headers:
						if h.find(inName)<0:
							length-=1
						else:
							break
					if length ==0:
						includes.add(line)
			elif line.startswith('typedef'):
				typedef+=line
			else:
				content+=line
		hfile.close()
		out.write('#include ' + header_include(name, os.path.dirname(header),std_include)+'\n')
	out.write('%}\n')
	out.write('%import <llvm/Pass.h>\n')
	out.write('%include <cpointer.i>\n')
	out.write('%include <carrays.i>\n')
	out.write('%include <stdint.i>\n')
	out.write('%include <stl.i>\n')
	out.write('%include <typemaps.i>\n')
	out.write('%include <std_set.i>\n')
	out.write('%include <std_map.i>\n')
	out.write('%include <std_string.i>\n')
	out.write('%include <std_pair.i>\n')
	out.write('%include <std_list.i>\n')
	out.write('%include <std_vector.i>\n')
	out.write('%include <std_deque.i>\n')
	out.write('%include <inttypes.i>\n')
	out.write('%include <attribute.i>\n')
	out.write('%include <exception.i>\n')
	out.write('%include <cdata.i>\n')
	out.write('%include <cmalloc.i>\n')
	out.write('%include <constraints.i>\n')
	out.write('%include <cwstring.i>\n')
	out.write('%include <intrusive_ptr.i>\n')
	out.write('%include <math.i>\n')
	out.write('%include <std_except.i>\n')
	out.write('%include <swigarch.i>\n')
	out.write('%include <swigrun.i>\n')
	out.write('%include <wchar.i>\n')
	out.write('%include <shared_ptr.i>\n')
	out.write('%pointer_class(PAG, PAGP)\n')
	out.write('%pointer_class(PointerAnalysis, PointerAnalysisP)\n')
	out.write('%pointer_class(Value, ValueP)\n')
	out.write('%pointer_class(SVFGSCC, SVFGSCCP)\n')
	out.write('%pointer_class(DDAClient, DDAClientP)\n')
	out.write('%pointer_class(SVFGEdge, SVFGEdgeP)\n')
	out.write('%pointer_class(SVFG, SVFGP)\n')
	out.write('%pointer_class(PAGNode,PAGNodeP)\n')
	out.write('%pointer_class(Function,FunctionP)\n')
	out.write('%pointer_class(StoreSVFGNode,StoreSVFGNodeP)\n')
	out.write('%pointer_class(MemObj,MemObjP)\n')
	out.write('%pointer_class(AddrSVFGNode,AddrSVFGNodeP)\n')
	out.write('%pointer_class(DDAStat,DDAStatP)\n')
	out.write('%pointer_class(Instruction,InstructionP)\n')
	out.write('%pointer_class(GepSVFGNode,GepSVFGNodeP)\n')
	out.write('%pointer_class(AndersenWaveDiff,AndersenWaveDiffP)\n')
	out.write('%pointer_class(PTACallGraph,PTACallGraphP)\n')
	out.write('%pointer_class(DirectSVFGEdge,DirectSVFGEdgeP)\n')
	out.write('%pointer_class(IndirectSVFGEdge,IndirectSVFGEdgeP)\n')
	out.write('%pointer_class(CallGraphSCC,CallGraphSCCP)\n')
	out.write('%pointer_class(LoadSVFGNode,LoadSVFGNodeP)\n')
	out.write('%pointer_class(AllocaInst,AllocaInstP)\n')
	out.write('%pointer_class(DiffPTDataTy,DiffPTDataTyP)\n')
	out.write('%pointer_class(PTAStat,PTAStatP)\n')
	out.write('%pointer_class(PTACallGraphNode,PTACallGraphNodeP)\n')
	out.write('%pointer_class(ICFG,ICFGP)\n')
	out.write('%pointer_class(PTDataTy,PTDataTyP)\n')
	out.write('%pointer_class(CallInst,CallInstP)\n')
	out.write('%pointer_class(IncDFPTDataTy,IncDFPTDataTyP)\n')
	out.write('%pointer_class(CHGraph,CHGraphP)\n')
	out.write('%pointer_class(TypeSystem,TypeSystemP)\n')
	out.write('%pointer_class(LLVMModuleSet,LLVMModuleSetP)\n')
	out.write('%pointer_class(LLVMContext,LLVMContextP)\n')
	out.write('%pointer_class(Module,ModuleP)\n')
	out.write('%pointer_class(GlobalAlias,GlobalAliasP)\n')
	out.write('%pointer_class(GlobalVariable,GlobalVariableP)\n')
	out.write('%pointer_class(u32_t,u32_tP)\n')
	out.write('%pointer_class(ConstantExpr,ConstantExprP)\n')
	out.write('%pointer_class(PointerType,PointerTypeP)\n')
	out.write('%pointer_class(BasicBlock,BasicBlockP)\n')
	out.write('%pointer_class(DominatorTree,DominatorTreeP)\n')

	out.write('namespace std {\n' +
			  '%template(DPTItemSet) set<DPIm>;\n' +
			  '%template(ConstSVFGEdgeSet) set<const SVFGEdge*>;\n' +
			  '%template(DPImToCPtSetMap) map<DPIm, CPtSet>;\n' +
			  '%template(DPMToCVarMap) map<DPIm,CVar>;\n' +
			  '%template(DPMToDPMMap) map<DPIm,DPIm>;\n' +
			  '%template(StoreToPMSetMap) map<const SVFGNode*, DPTItemSet>;\n' +
			  '%template(CallSiteSet) set<CallSite>;\n' +
			  '%template(FunctionSet) set<const Function*>;\n' +
			  '%template(VTableSet) set<const GlobalValue*>;\n' +
			  '%template(CallEdgeMap) map<CallSite, FunctionSet>;\n' +
			  '%template(PtrToBVPtsMap) map<NodeID,PointsTo>;\n' +
			  '%template(PtrCPtsMap) map<NodeID,CPtSet>;\n' +
			  '%template(NodePair) pair<NodeID, NodeID>;\n' +
			  '%template(NodeSet) set<NodeID>;\n' +
			  '%template(NodeVector) vector<NodeID>;\n' +
			  '%template(EdgeVector) vector<EdgeID>;\n' +
			  '%template(NodeList) list<NodeID>;\n' +
			  '%template(NodeDeque) deque<NodeID>;\n' +
			  '%template(NUMStatMap) map<const char*,u32_t>;\n' +
			  '%template(TIMEStatMap) map<const char*,double>;\n' +
			  '%template(FunctionSetType) vector<Function*> ;\n' +
			  '%template(GlobalSetType) vector<GlobalVariable*> ;\n' +
			  '%template(AliasSetType) vector<GlobalAlias*> ;\n' +
			  '%template(FunDeclToDefMapTy) map<const Function*, Function*>;\n' +
			  '%template(FunDefToDeclsMapTy) map<const Function*, FunctionSetType> ;\n' +
			  '%template(GlobalDefToRepMapTy) map<const GlobalVariable*, GlobalVariable*> ;\n' +
			  '%template(StringVector) vector<std::string>;\n' +
			  '%template(BasicBlockVector) vector<const BasicBlock*>;\n' +
			  '%template(InstructionVector) vector<const Instruction*>;\n' +
			  '}\n')
	for include in sys_includes:
		out.write(include)
	for include in includes:
		out.write(include)
	out.write(typedef)
	out.write(content)
	out.close()

if __name__ == '__main__':
	# if os.path.isdir("include") == True:
	# 	shutil.rmtree("include/")
	if os.path.isdir("src") == True:
		shutil.rmtree("src/")
	if os.path.isdir("jni") == True:
		shutil.rmtree("jni/")
	# generate_swig_interfaces(sys.argv[1])
	# print 'Finish interface files generation'
	# generate_swig_single(sys.argv[1])
	regenerate_swig_interfaces(sys.argv[2])
	print 'Generate wrap file'
