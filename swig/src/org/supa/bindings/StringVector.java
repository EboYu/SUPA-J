/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class StringVector extends java.util.AbstractList<String> implements java.util.RandomAccess {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected StringVector(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(StringVector obj) {
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
        SUPAJNI.delete_StringVector(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public StringVector(String[] initialElements) {
    this();
    reserve(initialElements.length);

    for (String element : initialElements) {
      add(element);
    }
  }

  public StringVector(Iterable<String> initialElements) {
    this();
    for (String element : initialElements) {
      add(element);
    }
  }

  public String get(int index) {
    return doGet(index);
  }

  public String set(int index, String e) {
    return doSet(index, e);
  }

  public boolean add(String e) {
    modCount++;
    doAdd(e);
    return true;
  }

  public void add(int index, String e) {
    modCount++;
    doAdd(index, e);
  }

  public String remove(int index) {
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

  public StringVector() {
    this(SUPAJNI.new_StringVector__SWIG_0(), true);
  }

  public StringVector(StringVector other) {
    this(SUPAJNI.new_StringVector__SWIG_1(StringVector.getCPtr(other), other), true);
  }

  public long capacity() {
    return SUPAJNI.StringVector_capacity(swigCPtr, this);
  }

  public void reserve(long n) {
    SUPAJNI.StringVector_reserve(swigCPtr, this, n);
  }

  public boolean isEmpty() {
    return SUPAJNI.StringVector_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.StringVector_clear(swigCPtr, this);
  }

  public StringVector(int count, String value) {
    this(SUPAJNI.new_StringVector__SWIG_2(count, value), true);
  }

  private int doSize() {
    return SUPAJNI.StringVector_doSize(swigCPtr, this);
  }

  private void doAdd(String x) {
    SUPAJNI.StringVector_doAdd__SWIG_0(swigCPtr, this, x);
  }

  private void doAdd(int index, String x) {
    SUPAJNI.StringVector_doAdd__SWIG_1(swigCPtr, this, index, x);
  }

  private String doRemove(int index) {
    return SUPAJNI.StringVector_doRemove(swigCPtr, this, index);
  }

  private String doGet(int index) {
    return SUPAJNI.StringVector_doGet(swigCPtr, this, index);
  }

  private String doSet(int index, String val) {
    return SUPAJNI.StringVector_doSet(swigCPtr, this, index, val);
  }

  private void doRemoveRange(int fromIndex, int toIndex) {
    SUPAJNI.StringVector_doRemoveRange(swigCPtr, this, fromIndex, toIndex);
  }

}