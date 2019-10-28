/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class NodeVector extends java.util.AbstractList<Long> implements java.util.RandomAccess {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected NodeVector(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(NodeVector obj) {
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
        SUPAJNI.delete_NodeVector(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public NodeVector(long[] initialElements) {
    this();
    reserve(initialElements.length);

    for (long element : initialElements) {
      add(element);
    }
  }

  public NodeVector(Iterable<Long> initialElements) {
    this();
    for (long element : initialElements) {
      add(element);
    }
  }

  public Long get(int index) {
    return doGet(index);
  }

  public Long set(int index, Long e) {
    return doSet(index, e);
  }

  public boolean add(Long e) {
    modCount++;
    doAdd(e);
    return true;
  }

  public void add(int index, Long e) {
    modCount++;
    doAdd(index, e);
  }

  public Long remove(int index) {
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

  public NodeVector() {
    this(SUPAJNI.new_NodeVector__SWIG_0(), true);
  }

  public NodeVector(NodeVector other) {
    this(SUPAJNI.new_NodeVector__SWIG_1(NodeVector.getCPtr(other), other), true);
  }

  public long capacity() {
    return SUPAJNI.NodeVector_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    SUPAJNI.NodeVector_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return SUPAJNI.NodeVector_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.NodeVector_clear(swigCPtr, this);
  }

  public NodeVector(int count, long value) {
    this(SUPAJNI.new_NodeVector__SWIG_2(count, value), true);
  }

  private int doSize() {
    return SUPAJNI.NodeVector_doSize(swigCPtr, this);
  }

  private void doAdd(long x) {
    SUPAJNI.NodeVector_doAdd__SWIG_0(swigCPtr, this, x);
  }

  private void doAdd(int index, long x) {
    SUPAJNI.NodeVector_doAdd__SWIG_1(swigCPtr, this, index, x);
  }

  private long doRemove(int index) {
    return SUPAJNI.NodeVector_doRemove(swigCPtr, this, index);
  }

  private long doGet(int index) {
    return SUPAJNI.NodeVector_doGet(swigCPtr, this, index);
  }

  private long doSet(int index, long val) {
    return SUPAJNI.NodeVector_doSet(swigCPtr, this, index, val);
  }

  private void doRemoveRange(int fromIndex, int toIndex) {
    SUPAJNI.NodeVector_doRemoveRange(swigCPtr, this, fromIndex, toIndex);
  }

}