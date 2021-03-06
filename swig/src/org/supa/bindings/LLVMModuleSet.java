/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class LLVMModuleSet {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected LLVMModuleSet(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(LLVMModuleSet obj) {
    return (obj == null) ? 0 : obj.swigCPtr;
  }

  @SuppressWarnings("deprecation")
  protected void finalize() {
    delete();
  }

  public synchronized void delete() {
    if (swigCPtr != 0) {
      if (swigCMemOwn) {
        swigCMemOwn = false;
        SUPAJNI.delete_LLVMModuleSet(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public LLVMModuleSet(StringVector moduleNameVec) {
    this(SUPAJNI.new_LLVMModuleSetMNV(StringVector.getCPtr(moduleNameVec), moduleNameVec), true);
  }

  public LLVMModuleSet(SWIGTYPE_p_llvm__Module mod) {
    this(SUPAJNI.new_LLVMModuleSetMOD(SWIGTYPE_p_llvm__Module.getCPtr(mod)), true);
  }

  public LLVMModuleSet(SWIGTYPE_p_llvm__Module mod) {
    this(SUPAJNI.new_LLVMModuleSetAMOD(SWIGTYPE_p_llvm__Module.getCPtr(mod)), true);
  }

  public LLVMModuleSet() {
    this(SUPAJNI.new_LLVMModuleSet(), true);
  }

  public void build(StringVector moduleNameVec) {
    SUPAJNI.LLVMModuleSet_build(swigCPtr, this, StringVector.getCPtr(moduleNameVec), moduleNameVec);
  }

  public long getModuleNum() {
    return SUPAJNI.LLVMModuleSet_getModuleNum(swigCPtr, this);
  }

  public SWIGTYPE_p_llvm__Module getModule(long idx) {
    long cPtr = SUPAJNI.LLVMModuleSet_getModule(swigCPtr, this, idx);
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Module(cPtr, false);
  }

  public SWIGTYPE_p_llvm__Module getModuleRef(long idx) {
    return new SWIGTYPE_p_llvm__Module(SUPAJNI.LLVMModuleSet_getModuleRef(swigCPtr, this, idx), false);
  }

  public void dumpModulesToFile(String suffix) {
    SUPAJNI.LLVMModuleSet_dumpModulesToFile(swigCPtr, this, suffix);
  }

  public boolean hasDefinition(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.LLVMModuleSet_hasDefinition(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public SWIGTYPE_p_llvm__Function getDefinition(SWIGTYPE_p_llvm__Function fun) {
    long cPtr = SUPAJNI.LLVMModuleSet_getDefinition(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(fun));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
  }

  public boolean hasDeclaration(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.LLVMModuleSet_hasDeclaration(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public FunctionSetType getDeclaration(SWIGTYPE_p_llvm__Function fun) {
    return new FunctionSetType(SUPAJNI.LLVMModuleSet_getDeclaration(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(fun)), false);
  }

  public boolean hasGlobalRep(SWIGTYPE_p_llvm__GlobalVariable val) {
    return SUPAJNI.LLVMModuleSet_hasGlobalRep(swigCPtr, this, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(val));
  }

  public SWIGTYPE_p_llvm__GlobalVariable getGlobalRep(SWIGTYPE_p_llvm__GlobalVariable val) {
    long cPtr = SUPAJNI.LLVMModuleSet_getGlobalRep(swigCPtr, this, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__GlobalVariable(cPtr, false);
  }

  public SWIGTYPE_p_std__vectorT_llvm__Function_p_t__iterator begin() {
    return new SWIGTYPE_p_std__vectorT_llvm__Function_p_t__iterator(SUPAJNI.LLVMModuleSet_begin(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__Function_p_t__const_iterator begin_const() {
    return new SWIGTYPE_p_std__vectorT_llvm__Function_p_t__const_iterator(SUPAJNI.LLVMModuleSet_begin_const(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__Function_p_t__iterator end() {
    return new SWIGTYPE_p_std__vectorT_llvm__Function_p_t__iterator(SUPAJNI.LLVMModuleSet_end(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__Function_p_t__const_iterator end_const() {
    return new SWIGTYPE_p_std__vectorT_llvm__Function_p_t__const_iterator(SUPAJNI.LLVMModuleSet_end_const(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__iterator global_begin() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__iterator(SUPAJNI.LLVMModuleSet_global_begin(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__const_iterator global_begin_const() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__const_iterator(SUPAJNI.LLVMModuleSet_global_begin_const(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__iterator global_end() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__iterator(SUPAJNI.LLVMModuleSet_global_end(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__const_iterator global_end_const() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalVariable_p_t__const_iterator(SUPAJNI.LLVMModuleSet_global_end_const(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__iterator alias_begin() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__iterator(SUPAJNI.LLVMModuleSet_alias_begin(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__const_iterator alias_begin_const() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__const_iterator(SUPAJNI.LLVMModuleSet_alias_begin_const(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__iterator alias_end() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__iterator(SUPAJNI.LLVMModuleSet_alias_end(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__const_iterator alias_end_const() {
    return new SWIGTYPE_p_std__vectorT_llvm__GlobalAlias_p_t__const_iterator(SUPAJNI.LLVMModuleSet_alias_end_const(swigCPtr, this), true);
  }

}
