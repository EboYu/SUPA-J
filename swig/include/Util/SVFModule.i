%module SVFModule_SWIG
%{
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/SVFModule.h"
%}
%include <cpointer.i>
%include <stdint.i>
%include <std_vector.i>
%include <std_map.i>
%include <std_shared_ptr.i>
%include <std_string.i>
%pointer_class(LLVMModuleSet,LLVMModuleSetP)
%pointer_class(Function,FunctionP)
%pointer_class(LLVMContext,LLVMContextP)
%pointer_class(Module,ModuleP)
%pointer_class(GlobalAlias,GlobalAliasP)
%pointer_class(GlobalVariable,GlobalVariableP)
namespace std {
%template(FunctionSetType) vector<Function*> ;
%template(GlobalSetType) vector<GlobalVariable*> ;
%template(AliasSetType) vector<GlobalAlias*> ;
%template(FunDeclToDefMapTy) map<const Function*, Function*>;
%template(FunDefToDeclsMapTy) map<const Function*, FunctionSetType> ;
%template(GlobalDefToRepMapTy) map<const GlobalVariable*, GlobalVariable*> ;
%template(StringVector) vector<std::string>;
}
%copyctor LLVMModuleSet;
//===- SVFModule.h -- SVFModule class-----------------------------------------//
//
//                     SVF: Static Value-Flow Analysis
//
// Copyright (C) <2013-2017>  <Yulei Sui>
//

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//===----------------------------------------------------------------------===//

/*
 * SVFModule.h
 *
 *  Created on: Aug 4, 2017
 *      Author: Xiaokang Fan
 */

#ifndef SVFMODULE_H_
#define SVFMODULE_H_

#include "Util/BasicTypes.h"

class LLVMModuleSet {
public:
    typedef std::vector<Function*> FunctionSetType;
    typedef std::vector<GlobalVariable*> GlobalSetType;
    typedef std::vector<GlobalAlias*> AliasSetType;

    typedef std::map<const Function*, Function*> FunDeclToDefMapTy;
    typedef std::map<const Function*, FunctionSetType> FunDefToDeclsMapTy;
    typedef std::map<const GlobalVariable*, GlobalVariable*> GlobalDefToRepMapTy;

    /// Iterators type def
    typedef FunctionSetType::iterator iterator;
    typedef FunctionSetType::const_iterator const_iterator;
    typedef GlobalSetType::iterator global_iterator;
    typedef GlobalSetType::const_iterator const_global_iterator;
    typedef AliasSetType::iterator alias_iterator;
    typedef AliasSetType::const_iterator const_alias_iterator;

private:
    u32_t moduleNum;
    LLVMContext *cxts;
    std::unique_ptr<Module> *modules;

    FunctionSetType FunctionSet;  ///< The Functions in the module
    GlobalSetType GlobalSet;      ///< The Global Variables in the module
    AliasSetType AliasSet;        ///< The Aliases in the module

    /// Function declaration to function definition map
    FunDeclToDefMapTy FunDeclToDefMap;
    /// Function definition to function declaration map
    FunDefToDeclsMapTy FunDefToDeclsMap;
    /// Global definition to a rep definition map
    GlobalDefToRepMapTy GlobalDefToRepMap;

public:
    /// Constructor
    %rename(LLVMModuleSetMOD) LLVMModuleSet(Module *mod);
    %rename(LLVMModuleSetAMOD) LLVMModuleSet(Module &mod);
    %rename(LLVMModuleSetMNV) LLVMModuleSet(const std::vector<std::string> &moduleNameVec);
    LLVMModuleSet() {}
    
    LLVMModuleSet(const std::vector<std::string> &moduleNameVec);
    LLVMModuleSet(Module *mod);
    LLVMModuleSet(Module &mod);
    

    void build(const std::vector<std::string> &moduleNameVec);

    u32_t getModuleNum() const {
        return moduleNum;
    }

    Module *getModule(u32_t idx) const {
        assert(idx < moduleNum && "Out of range.");
        return modules[idx].get();
    }

    Module &getModuleRef(u32_t idx) const {
        assert(idx < moduleNum && "Out of range.");
        return *(modules[idx].get());
    }

    // Dump modules to files
    void dumpModulesToFile(const std::string suffix);

    /// Fun decl --> def
    bool hasDefinition(const Function *fun) const {
        assert(fun->isDeclaration() && "not a function declaration?");
        FunDeclToDefMapTy::const_iterator it = FunDeclToDefMap.find(fun);
        return it != FunDeclToDefMap.end();
    }

    Function *getDefinition(const Function *fun) const {
        assert(fun->isDeclaration() && "not a function declaration?");
        FunDeclToDefMapTy::const_iterator it = FunDeclToDefMap.find(fun);
        assert(it != FunDeclToDefMap.end() && "has no definition?");
        return it->second;
    }

    /// Fun def --> decl
    bool hasDeclaration(const Function *fun) const {
        assert(!fun->isDeclaration() && "not a function definition?");
        FunDefToDeclsMapTy::const_iterator it = FunDefToDeclsMap.find(fun);
        return it != FunDefToDeclsMap.end();
    }

    const FunctionSetType &getDeclaration(const Function *fun) const {
        assert(!fun->isDeclaration() && "not a function definition?");
        FunDefToDeclsMapTy::const_iterator it = FunDefToDeclsMap.find(fun);
        assert(it != FunDefToDeclsMap.end() && "has no declaration?");
        return it->second;
    }

    /// Global to rep
    bool hasGlobalRep(const GlobalVariable *val) const {
        GlobalDefToRepMapTy::const_iterator it = GlobalDefToRepMap.find(val);
        return it != GlobalDefToRepMap.end();
    }

    GlobalVariable *getGlobalRep(const GlobalVariable *val) const {
        GlobalDefToRepMapTy::const_iterator it = GlobalDefToRepMap.find(val);
        assert(it != GlobalDefToRepMap.end() && "has no rep?");
        return it->second;
    }

    /// Iterators
    ///@{
    iterator begin() {
        return FunctionSet.begin();
    }
    %rename(begin_const) begin() const;
    const_iterator begin() const {
        return FunctionSet.begin();
    }
    iterator end() {
        return FunctionSet.end();
    }
    %rename(end_const) end() const;
    const_iterator end() const {
        return FunctionSet.end();
    }

    global_iterator global_begin() {
        return GlobalSet.begin();
    }
    %rename(global_begin_const) global_begin() const;
    const_global_iterator global_begin() const {
        return GlobalSet.begin();
    }
    global_iterator global_end() {
        return GlobalSet.end();
    }
    %rename(global_end_const) global_end() const;
    const_global_iterator global_end() const {
        return GlobalSet.end();
    }

    alias_iterator alias_begin() {
        return AliasSet.begin();
    }

%rename(alias_begin_const) alias_begin() const;
    const_alias_iterator alias_begin() const {
        return AliasSet.begin();
    }
    alias_iterator alias_end() {
        return AliasSet.end();
    }
%rename(alias_end_const) alias_end() const;
    const_alias_iterator alias_end() const {
        return AliasSet.end();
    }
    ///@}

private:
    void loadModules(const std::vector<std::string> &moduleNameVec);
    void addSVFMain();
    void initialize();
    void buildFunToFunMap();
    void buildGlobalDefToRepMap();
};

class SVFModule {
public:
    typedef LLVMModuleSet::FunctionSetType FunctionSetType;
    typedef LLVMModuleSet::GlobalSetType GlobalSetType;
    typedef LLVMModuleSet::AliasSetType AliasSetType;

    typedef LLVMModuleSet::FunDeclToDefMapTy FunDeclToDefMapTy;
    typedef LLVMModuleSet::FunDefToDeclsMapTy FunDefToDeclsMapTy;
    typedef LLVMModuleSet::GlobalDefToRepMapTy GlobalDefToRepMapTy;

    /// Iterators type def
    typedef FunctionSetType::iterator iterator;
    typedef FunctionSetType::const_iterator const_iterator;
    typedef GlobalSetType::iterator global_iterator;
    typedef GlobalSetType::const_iterator const_global_iterator;
    typedef AliasSetType::iterator alias_iterator;
    typedef AliasSetType::const_iterator const_alias_iterator;

private:
    static LLVMModuleSet *llvmModuleSet;
    static std::string pagReadFromTxt;

public:
    /// Constructors
    %rename(SVFModulePointer) SVFModule(Module *mod);
    %rename(SVFModuleAddr) SVFModule(Module &mod);
    %rename(SVFModuleMNV) SVFModule(const std::vector<std::string> &moduleNameVec);
    
    
    SVFModule(const std::vector<std::string> &moduleNameVec) {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet(moduleNameVec);
    }
    SVFModule(Module *mod) {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet(mod);
    }
    SVFModule(Module &mod) {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet(mod);
    }
    SVFModule() {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet;
    }

    

    static inline LLVMModuleSet *getLLVMModuleSet() {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet;
        return llvmModuleSet;
    }

    static inline void setPagFromTXT(std::string txt) {
        pagReadFromTxt = txt;
    }

    static inline std::string pagFileName() {
        return pagReadFromTxt;
    }

    static inline bool pagReadFromTXT() {
    		if(pagReadFromTxt.empty())
    			return false;
    		else
    			return true;
    }

    static void releaseLLVMModuleSet() {
        if (llvmModuleSet)
            delete llvmModuleSet;
        llvmModuleSet = NULL;
    }

    bool empty() const {
        return getModuleNum() == 0;
    }

    /// Methods from LLVMModuleSet
    u32_t getModuleNum() const {
        return llvmModuleSet->getModuleNum();
    }

    Module *getModule(u32_t idx) const {
        return llvmModuleSet->getModule(idx);
    }

    Module &getModuleRef(u32_t idx) const {
        return llvmModuleSet->getModuleRef(idx);
    }

    // Dump modules to files
    void dumpModulesToFile(const std::string suffix) const {
        llvmModuleSet->dumpModulesToFile(suffix);
    }

    /// Fun decl --> def
    bool hasDefinition(const Function *fun) const {
        return llvmModuleSet->hasDefinition(fun);
    }

    Function *getDefinition(const Function *fun) const {
        return llvmModuleSet->getDefinition(fun);
    }

    /// Fun def --> decl
    bool hasDeclaration(const Function *fun) const {
        return llvmModuleSet->hasDeclaration(fun);
    }

    const FunctionSetType &getDeclaration(const Function *fun) const {
        return llvmModuleSet->getDeclaration(fun);
    }

    /// Global to rep
    bool hasGlobalRep(const GlobalVariable *val) const {
        return llvmModuleSet->hasGlobalRep(val);
    }

    GlobalVariable *getGlobalRep(const GlobalVariable *val) const {
        return llvmModuleSet->getGlobalRep(val);
    }

    /// Iterators
    ///@{
    iterator begin() {
        return llvmModuleSet->begin();
    }
    %rename(begin_const) begin() const;
    const_iterator begin() const {
        return llvmModuleSet->begin();
    }
    iterator end() {
        return llvmModuleSet->end();
    }
    %rename(end_const) end() const;
    const_iterator end() const {
        return llvmModuleSet->end();
    }

    Module *getMainLLVMModule() const {
        return llvmModuleSet->getModule(0);
    }

	const std::string& getModuleIdentifier() const {
		if (pagReadFromTxt.empty()) {
			assert(!empty() && "empty LLVM module!!");
			return getMainLLVMModule()->getModuleIdentifier();
		} else {
			return pagReadFromTxt;
		}
	}

    LLVMContext& getContext() const {
        assert(!empty() && "empty LLVM module!!");
        return getMainLLVMModule()->getContext();
    }

    inline Function* getFunction(StringRef name) const {
        Function* fun = NULL;
        for (u32_t i = 0; i < getModuleNum(); ++i) {
            Module *mod = llvmModuleSet->getModule(i);
            fun = mod->getFunction(name);
            if(fun && !fun->isDeclaration()) {
                return fun;
            }
        }
        return fun;
    }

%rename(global_begin_const) global_begin() const;
%rename(global_end_const) global_end() const;
%rename(alias_begin_const) alias_begin() const;
%rename(alias_end_const) alias_end() const;
    global_iterator global_begin() {
        return llvmModuleSet->global_begin();
    }
    const_global_iterator global_begin() const {
        return llvmModuleSet->global_begin();
    }
    global_iterator global_end() {
        return llvmModuleSet->global_end();
    }
    const_global_iterator global_end() const {
        return llvmModuleSet->global_end();
    }

    alias_iterator alias_begin() {
        return llvmModuleSet->alias_begin();
    }
    const_alias_iterator alias_begin() const {
        return llvmModuleSet->alias_begin();
    }
    alias_iterator alias_end() {
        return llvmModuleSet->alias_end();
    }
    const_alias_iterator alias_end() const {
        return llvmModuleSet->alias_end();
    }
    ///@}
};


#endif /* SVFMODULE_H_ */
