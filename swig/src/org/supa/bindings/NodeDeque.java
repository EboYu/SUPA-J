/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class NodeDeque {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected NodeDeque(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(NodeDeque obj) {
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
        SUPAJNI.delete_NodeDeque(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public boolean empty() {
    return SUPAJNI.NodeDeque_empty(swigCPtr, this);
  }

  public NodeDeque() {
    this(SUPAJNI.new_NodeDeque__SWIG_0(), true);
  }

  public NodeDeque(long size, long value) {
    this(SUPAJNI.new_NodeDeque__SWIG_1(size, value), true);
  }

  public NodeDeque(long size) {
    this(SUPAJNI.new_NodeDeque__SWIG_2(size), true);
  }

  public NodeDeque(NodeDeque arg0) {
    this(SUPAJNI.new_NodeDeque__SWIG_3(NodeDeque.getCPtr(arg0), arg0), true);
  }

  public void assign(long n, long value) {
    SUPAJNI.NodeDeque_assign(swigCPtr, this, n, value);
  }

  public void swap(NodeDeque x) {
    SUPAJNI.NodeDeque_swap(swigCPtr, this, NodeDeque.getCPtr(x), x);
  }

  public long size() {
    return SUPAJNI.NodeDeque_size(swigCPtr, this);
  }

  public long max_size() {
    return SUPAJNI.NodeDeque_max_size(swigCPtr, this);
  }

  public void resize(long n, long c) {
    SUPAJNI.NodeDeque_resize__SWIG_0(swigCPtr, this, n, c);
  }

  public void resize(long n) {
    SUPAJNI.NodeDeque_resize__SWIG_1(swigCPtr, this, n);
  }

  public long front() {
    return SUPAJNI.NodeDeque_front(swigCPtr, this);
  }

  public long back() {
    return SUPAJNI.NodeDeque_back(swigCPtr, this);
  }

  public void push_front(long x) {
    SUPAJNI.NodeDeque_push_front(swigCPtr, this, x);
  }

  public void push_back(long x) {
    SUPAJNI.NodeDeque_push_back(swigCPtr, this, x);
  }

  public void pop_front() {
    SUPAJNI.NodeDeque_pop_front(swigCPtr, this);
  }

  public void pop_back() {
    SUPAJNI.NodeDeque_pop_back(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.NodeDeque_clear(swigCPtr, this);
  }

  public long getitem(int i) {
    return SUPAJNI.NodeDeque_getitem(swigCPtr, this, i);
  }

  public void setitem(int i, long x) {
    SUPAJNI.NodeDeque_setitem(swigCPtr, this, i, x);
  }

  public void delitem(int i) {
    SUPAJNI.NodeDeque_delitem(swigCPtr, this, i);
  }

  public NodeDeque getslice(int i, int j) {
    return new NodeDeque(SUPAJNI.NodeDeque_getslice(swigCPtr, this, i, j), true);
  }

  public void setslice(int i, int j, NodeDeque v) {
    SUPAJNI.NodeDeque_setslice(swigCPtr, this, i, j, NodeDeque.getCPtr(v), v);
  }

  public void delslice(int i, int j) {
    SUPAJNI.NodeDeque_delslice(swigCPtr, this, i, j);
  }

}