/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class GlobalSetType extends java.util.AbstractList<SWIGTYPE_p_llvm__GlobalVariable> implements java.util.RandomAccess {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected GlobalSetType(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(GlobalSetType obj) {
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
        SUPAJNI.delete_GlobalSetType(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public GlobalSetType(SWIGTYPE_p_llvm__GlobalVariable[] initialElements) {
    this();
    reserve(initialElements.length);

    for (SWIGTYPE_p_llvm__GlobalVariable element : initialElements) {
      add(element);
    }
  }

  public GlobalSetType(Iterable<SWIGTYPE_p_llvm__GlobalVariable> initialElements) {
    this();
    for (SWIGTYPE_p_llvm__GlobalVariable element : initialElements) {
      add(element);
    }
  }

  public SWIGTYPE_p_llvm__GlobalVariable get(int index) {
    return doGet(index);
  }

  public SWIGTYPE_p_llvm__GlobalVariable set(int index, SWIGTYPE_p_llvm__GlobalVariable e) {
    return doSet(index, e);
  }

  public boolean add(SWIGTYPE_p_llvm__GlobalVariable e) {
    modCount++;
    doAdd(e);
    return true;
  }

  public void add(int index, SWIGTYPE_p_llvm__GlobalVariable e) {
    modCount++;
    doAdd(index, e);
  }

  public SWIGTYPE_p_llvm__GlobalVariable remove(int index) {
    modCount++;
    return doRemove(index);
  }

  protected void removeRange(int fromIndex, int toIndex) {
    modCount++;
    doRemoveRange(fromIndex, toIndex);
  }

  public int size() {
    return doSize();
  }

  public GlobalSetType() {
    this(SUPAJNI.new_GlobalSetType__SWIG_0(), true);
  }

  public GlobalSetType(GlobalSetType other) {
    this(SUPAJNI.new_GlobalSetType__SWIG_1(GlobalSetType.getCPtr(other), other), true);
  }

  public long capacity() {
    return SUPAJNI.GlobalSetType_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    SUPAJNI.GlobalSetType_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return SUPAJNI.GlobalSetType_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.GlobalSetType_clear(swigCPtr, this);
  }

  public GlobalSetType(int count, SWIGTYPE_p_llvm__GlobalVariable value) {
    this(SUPAJNI.new_GlobalSetType__SWIG_2(count, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(value)), true);
  }

  private int doSize() {
    return SUPAJNI.GlobalSetType_doSize(swigCPtr, this);
  }

  private void doAdd(SWIGTYPE_p_llvm__GlobalVariable x) {
    SUPAJNI.GlobalSetType_doAdd__SWIG_0(swigCPtr, this, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(x));
  }

  private void doAdd(int index, SWIGTYPE_p_llvm__GlobalVariable x) {
    SUPAJNI.GlobalSetType_doAdd__SWIG_1(swigCPtr, this, index, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(x));
  }

  private SWIGTYPE_p_llvm__GlobalVariable doRemove(int index) {
    long cPtr = SUPAJNI.GlobalSetType_doRemove(swigCPtr, this, index);
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__GlobalVariable(cPtr, false);
  }

  private SWIGTYPE_p_llvm__GlobalVariable doGet(int index) {
    long cPtr = SUPAJNI.GlobalSetType_doGet(swigCPtr, this, index);
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__GlobalVariable(cPtr, false);
  }

  private SWIGTYPE_p_llvm__GlobalVariable doSet(int index, SWIGTYPE_p_llvm__GlobalVariable val) {
    long cPtr = SUPAJNI.GlobalSetType_doSet(swigCPtr, this, index, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__GlobalVariable(cPtr, false);
  }

  private void doRemoveRange(int fromIndex, int toIndex) {
    SUPAJNI.GlobalSetType_doRemoveRange(swigCPtr, this, fromIndex, toIndex);
  }

}
