/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class InstructionVector extends java.util.AbstractList<SWIGTYPE_p_llvm__Instruction> implements java.util.RandomAccess {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected InstructionVector(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(InstructionVector obj) {
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
        SUPAJNI.delete_InstructionVector(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public InstructionVector(SWIGTYPE_p_llvm__Instruction[] initialElements) {
    this();
    reserve(initialElements.length);

    for (SWIGTYPE_p_llvm__Instruction element : initialElements) {
      add(element);
    }
  }

  public InstructionVector(Iterable<SWIGTYPE_p_llvm__Instruction> initialElements) {
    this();
    for (SWIGTYPE_p_llvm__Instruction element : initialElements) {
      add(element);
    }
  }

  public SWIGTYPE_p_llvm__Instruction get(int index) {
    return doGet(index);
  }

  public SWIGTYPE_p_llvm__Instruction set(int index, SWIGTYPE_p_llvm__Instruction e) {
    return doSet(index, e);
  }

  public boolean add(SWIGTYPE_p_llvm__Instruction e) {
    modCount++;
    doAdd(e);
    return true;
  }

  public void add(int index, SWIGTYPE_p_llvm__Instruction e) {
    modCount++;
    doAdd(index, e);
  }

  public SWIGTYPE_p_llvm__Instruction remove(int index) {
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

  public InstructionVector() {
    this(SUPAJNI.new_InstructionVector__SWIG_0(), true);
  }

  public InstructionVector(InstructionVector other) {
    this(SUPAJNI.new_InstructionVector__SWIG_1(InstructionVector.getCPtr(other), other), true);
  }

  public long capacity() {
    return SUPAJNI.InstructionVector_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    SUPAJNI.InstructionVector_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return SUPAJNI.InstructionVector_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.InstructionVector_clear(swigCPtr, this);
  }

  public InstructionVector(int count, SWIGTYPE_p_llvm__Instruction value) {
    this(SUPAJNI.new_InstructionVector__SWIG_2(count, SWIGTYPE_p_llvm__Instruction.getCPtr(value)), true);
  }

  private int doSize() {
    return SUPAJNI.InstructionVector_doSize(swigCPtr, this);
  }

  private void doAdd(SWIGTYPE_p_llvm__Instruction x) {
    SUPAJNI.InstructionVector_doAdd__SWIG_0(swigCPtr, this, SWIGTYPE_p_llvm__Instruction.getCPtr(x));
  }

  private void doAdd(int index, SWIGTYPE_p_llvm__Instruction x) {
    SUPAJNI.InstructionVector_doAdd__SWIG_1(swigCPtr, this, index, SWIGTYPE_p_llvm__Instruction.getCPtr(x));
  }

  private SWIGTYPE_p_llvm__Instruction doRemove(int index) {
    long cPtr = SUPAJNI.InstructionVector_doRemove(swigCPtr, this, index);
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Instruction(cPtr, false);
  }

  private SWIGTYPE_p_llvm__Instruction doGet(int index) {
    long cPtr = SUPAJNI.InstructionVector_doGet(swigCPtr, this, index);
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Instruction(cPtr, false);
  }

  private SWIGTYPE_p_llvm__Instruction doSet(int index, SWIGTYPE_p_llvm__Instruction val) {
    long cPtr = SUPAJNI.InstructionVector_doSet(swigCPtr, this, index, SWIGTYPE_p_llvm__Instruction.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Instruction(cPtr, false);
  }

  private void doRemoveRange(int fromIndex, int toIndex) {
    SUPAJNI.InstructionVector_doRemoveRange(swigCPtr, this, fromIndex, toIndex);
  }

}
