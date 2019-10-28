/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class ConstSVFGEdgeSet extends java.util.AbstractSet<SWIGTYPE_p_SVFGEdge> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected ConstSVFGEdgeSet(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(ConstSVFGEdgeSet obj) {
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
        SUPAJNI.delete_ConstSVFGEdgeSet(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public ConstSVFGEdgeSet(java.util.Collection<? extends SWIGTYPE_p_SVFGEdge> collection) {
    this();
    addAll(collection);
  }

  public int size() {
    return sizeImpl();
  }

  public boolean add(SWIGTYPE_p_SVFGEdge key) {
    return addImpl(key);
  }

  public boolean addAll(java.util.Collection<? extends SWIGTYPE_p_SVFGEdge> collection) {
    boolean didAddElement = false;
    for (java.lang.Object object : collection) {
      didAddElement |= add((SWIGTYPE_p_SVFGEdge)object);
    }

    return didAddElement;
  }

  public java.util.Iterator<SWIGTYPE_p_SVFGEdge> iterator() {
    return new java.util.Iterator<SWIGTYPE_p_SVFGEdge>() {
      private Iterator curr;
      private Iterator end;

      private java.util.Iterator<SWIGTYPE_p_SVFGEdge> init() {
        curr = ConstSVFGEdgeSet.this.begin();
        end = ConstSVFGEdgeSet.this.end();
        return this;
      }

      public SWIGTYPE_p_SVFGEdge next() {
        if (!hasNext()) {
          throw new java.util.NoSuchElementException();
        }

        // Save the current position, increment it,
        // then return the value at the position before the increment.
        final SWIGTYPE_p_SVFGEdge currValue = curr.derefUnchecked();
        curr.incrementUnchecked();
        return currValue;
      }

      public boolean hasNext() {
        return curr.isNot(end);
      }
    }.init();
  }

  public boolean containsAll(java.util.Collection<?> collection) {
    for (java.lang.Object object : collection) {
      if (!contains(object)) {
        return false;
      }
    }

    return true;
  }

  public boolean contains(java.lang.Object object) {
    if (!(object instanceof SWIGTYPE_p_SVFGEdge)) {
      return false;
    }

    return containsImpl((SWIGTYPE_p_SVFGEdge)object);
  }

  public boolean removeAll(java.util.Collection<?> collection) {
    boolean didRemoveElement = false;
    for (java.lang.Object object : collection) {
      didRemoveElement |= remove(object);
    }

    return didRemoveElement;
  }

  public boolean remove(java.lang.Object object) {
    if (!(object instanceof SWIGTYPE_p_SVFGEdge)) {
      return false;
    }

    return removeImpl((SWIGTYPE_p_SVFGEdge)object);
  }

  static protected class Iterator {
    private transient long swigCPtr;
    protected transient boolean swigCMemOwn;
  
    protected Iterator(long cPtr, boolean cMemoryOwn) {
      swigCMemOwn = cMemoryOwn;
      swigCPtr = cPtr;
    }
  
    protected static long getCPtr(Iterator obj) {
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
          SUPAJNI.delete_ConstSVFGEdgeSet_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private void incrementUnchecked() {
      SUPAJNI.ConstSVFGEdgeSet_Iterator_incrementUnchecked(swigCPtr, this);
    }
  
    private SWIGTYPE_p_SVFGEdge derefUnchecked() {
      long cPtr = SUPAJNI.ConstSVFGEdgeSet_Iterator_derefUnchecked(swigCPtr, this);
      return (cPtr == 0) ? null : new SWIGTYPE_p_SVFGEdge(cPtr, false);
    }
  
    private boolean isNot(ConstSVFGEdgeSet.Iterator other) {
      return SUPAJNI.ConstSVFGEdgeSet_Iterator_isNot(swigCPtr, this, ConstSVFGEdgeSet.Iterator.getCPtr(other), other);
    }
  
  }

  public ConstSVFGEdgeSet() {
    this(SUPAJNI.new_ConstSVFGEdgeSet__SWIG_0(), true);
  }

  public ConstSVFGEdgeSet(ConstSVFGEdgeSet other) {
    this(SUPAJNI.new_ConstSVFGEdgeSet__SWIG_1(ConstSVFGEdgeSet.getCPtr(other), other), true);
  }

  public boolean isEmpty() {
    return SUPAJNI.ConstSVFGEdgeSet_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.ConstSVFGEdgeSet_clear(swigCPtr, this);
  }

  private ConstSVFGEdgeSet.Iterator begin() {
    return new ConstSVFGEdgeSet.Iterator(SUPAJNI.ConstSVFGEdgeSet_begin(swigCPtr, this), true);
  }

  private ConstSVFGEdgeSet.Iterator end() {
    return new ConstSVFGEdgeSet.Iterator(SUPAJNI.ConstSVFGEdgeSet_end(swigCPtr, this), true);
  }

  public boolean addImpl(SWIGTYPE_p_SVFGEdge key) {
    return SUPAJNI.ConstSVFGEdgeSet_addImpl(swigCPtr, this, SWIGTYPE_p_SVFGEdge.getCPtr(key));
  }

  private boolean containsImpl(SWIGTYPE_p_SVFGEdge key) {
    return SUPAJNI.ConstSVFGEdgeSet_containsImpl(swigCPtr, this, SWIGTYPE_p_SVFGEdge.getCPtr(key));
  }

  private boolean removeImpl(SWIGTYPE_p_SVFGEdge key) {
    return SUPAJNI.ConstSVFGEdgeSet_removeImpl(swigCPtr, this, SWIGTYPE_p_SVFGEdge.getCPtr(key));
  }

  private int sizeImpl() {
    return SUPAJNI.ConstSVFGEdgeSet_sizeImpl(swigCPtr, this);
  }

  private boolean hasNextImpl(ConstSVFGEdgeSet.Iterator itr) {
    return SUPAJNI.ConstSVFGEdgeSet_hasNextImpl(swigCPtr, this, ConstSVFGEdgeSet.Iterator.getCPtr(itr), itr);
  }

}
