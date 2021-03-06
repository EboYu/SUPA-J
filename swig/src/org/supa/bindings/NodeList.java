/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class NodeList extends java.util.AbstractSequentialList<Long> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected NodeList(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(NodeList obj) {
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
        SUPAJNI.delete_NodeList(swigCPtr);
      }
      swigCPtr = 0;
    }
  }

  public NodeList(java.util.Collection c) {
    this();
    java.util.ListIterator<Long> it = listIterator(0);
    // Special case the "copy constructor" here to avoid lots of cross-language calls
    for (java.lang.Object o : c) {
      it.add((Long)o);
    }
  }

  public int size() {
    return doSize();
  }

  public boolean add(Long value) {
    addLast(value);
    return true;
  }

  public java.util.ListIterator<Long> listIterator(int index) {
    return new java.util.ListIterator<Long>() {
      private Iterator pos;
      private Iterator last;

      private java.util.ListIterator<Long> init(int index) {
        if (index < 0 || index > NodeList.this.size())
          throw new IndexOutOfBoundsException("Index: " + index);
        pos = NodeList.this.begin();
	pos = pos.advance_unchecked(index);
        return this;
      }

      public void add(Long v) {
        // Technically we can invalidate last here, but this makes more sense
        last = NodeList.this.insert(pos, v);
      }

      public void set(Long v) {
        if (null == last) {
          throw new IllegalStateException();
        }
        last.set_unchecked(v);
      }

      public void remove() {
        if (null == last) {
          throw new IllegalStateException();
        }
        NodeList.this.remove(last);
        last = null;
      }

      public int previousIndex() {
        return NodeList.this.doPreviousIndex(pos);
      }

      public int nextIndex() {
        return NodeList.this.doNextIndex(pos);
      }

      public Long previous() {
        if (previousIndex() < 0) {
          throw new java.util.NoSuchElementException();
        }
        last = pos;
        pos = pos.previous_unchecked();
        return last.deref_unchecked();
      }

      public Long next() {
        if (!hasNext()) {
          throw new java.util.NoSuchElementException();
        }
        last = pos;
        pos = pos.next_unchecked();
        return last.deref_unchecked();
      }

      public boolean hasPrevious() {
        // This call to previousIndex() will be much slower than the hasNext() implementation, but it's simpler like this with C++ forward iterators
        return previousIndex() != -1;
      }

      public boolean hasNext() {
        return NodeList.this.doHasNext(pos);
      }
    }.init(index);
  }

  static public class Iterator {
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
          SUPAJNI.delete_NodeList_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    public void set_unchecked(long v) {
      SUPAJNI.NodeList_Iterator_set_unchecked(swigCPtr, v);
    }
  
    public SWIGTYPE_p_std__listT_unsigned_int_t__iterator next_unchecked() {
      return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_Iterator_next_unchecked(swigCPtr), true);
    }
  
    public SWIGTYPE_p_std__listT_unsigned_int_t__iterator previous_unchecked() {
      return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_Iterator_previous_unchecked(swigCPtr), true);
    }
  
    public long deref_unchecked() {
      return SUPAJNI.NodeList_Iterator_deref_unchecked(swigCPtr);
    }
  
    public SWIGTYPE_p_std__listT_unsigned_int_t__iterator advance_unchecked(long index) {
      return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_Iterator_advance_unchecked(swigCPtr, index), true);
    }
  
  }

  public NodeList() {
    this(SUPAJNI.new_NodeList__SWIG_0(), true);
  }

  public NodeList(NodeList other) {
    this(SUPAJNI.new_NodeList__SWIG_1(NodeList.getCPtr(other), other), true);
  }

  public boolean isEmpty() {
    return SUPAJNI.NodeList_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.NodeList_clear(swigCPtr, this);
  }

  public SWIGTYPE_p_std__listT_unsigned_int_t__iterator remove(SWIGTYPE_p_std__listT_unsigned_int_t__iterator pos) {
    return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_remove(swigCPtr, this, SWIGTYPE_p_std__listT_unsigned_int_t__iterator.getCPtr(pos)), true);
  }

  public void removeLast() {
    SUPAJNI.NodeList_removeLast(swigCPtr, this);
  }

  public void removeFirst() {
    SUPAJNI.NodeList_removeFirst(swigCPtr, this);
  }

  public void addLast(long value) {
    SUPAJNI.NodeList_addLast(swigCPtr, this, value);
  }

  public void addFirst(long value) {
    SUPAJNI.NodeList_addFirst(swigCPtr, this, value);
  }

  private SWIGTYPE_p_std__listT_unsigned_int_t__iterator begin() {
    return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_begin(swigCPtr, this), true);
  }

  public SWIGTYPE_p_std__listT_unsigned_int_t__iterator end() {
    return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_end(swigCPtr, this), true);
  }

  private SWIGTYPE_p_std__listT_unsigned_int_t__iterator insert(SWIGTYPE_p_std__listT_unsigned_int_t__iterator pos, long value) {
    return new SWIGTYPE_p_std__listT_unsigned_int_t__iterator(SUPAJNI.NodeList_insert(swigCPtr, this, SWIGTYPE_p_std__listT_unsigned_int_t__iterator.getCPtr(pos), value), true);
  }

  public NodeList(int count, long value) {
    this(SUPAJNI.new_NodeList__SWIG_2(count, value), true);
  }

  private int doSize() {
    return SUPAJNI.NodeList_doSize(swigCPtr, this);
  }

  private int doPreviousIndex(SWIGTYPE_p_std__listT_unsigned_int_t__iterator pos) {
    return SUPAJNI.NodeList_doPreviousIndex(swigCPtr, this, SWIGTYPE_p_std__listT_unsigned_int_t__iterator.getCPtr(pos));
  }

  private int doNextIndex(SWIGTYPE_p_std__listT_unsigned_int_t__iterator pos) {
    return SUPAJNI.NodeList_doNextIndex(swigCPtr, this, SWIGTYPE_p_std__listT_unsigned_int_t__iterator.getCPtr(pos));
  }

  private boolean doHasNext(SWIGTYPE_p_std__listT_unsigned_int_t__iterator pos) {
    return SUPAJNI.NodeList_doHasNext(swigCPtr, this, SWIGTYPE_p_std__listT_unsigned_int_t__iterator.getCPtr(pos));
  }

}
