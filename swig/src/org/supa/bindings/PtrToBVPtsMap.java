/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class PtrToBVPtsMap extends java.util.AbstractMap<Long, SWIGTYPE_p_llvm__SparseBitVectorT_t> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected PtrToBVPtsMap(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(PtrToBVPtsMap obj) {
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
        SUPAJNI.delete_PtrToBVPtsMap(swigCPtr);
      }
      swigCPtr = 0;
    }
  }


  public int size() {
    return sizeImpl();
  }

  public boolean containsKey(java.lang.Object key) {
    if (!(key instanceof Long)) {
      return false;
    }

    return containsImpl((Long)key);
  }

  public SWIGTYPE_p_llvm__SparseBitVectorT_t get(java.lang.Object key) {
    if (!(key instanceof Long)) {
      return null;
    }

    Iterator itr = find((Long) key);
    if (itr.isNot(end())) {
      return itr.getValue();
    }

    return null;
  }

  public SWIGTYPE_p_llvm__SparseBitVectorT_t put(Long key, SWIGTYPE_p_llvm__SparseBitVectorT_t value) {
    Iterator itr = find((Long) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_llvm__SparseBitVectorT_t oldValue = itr.getValue();
      itr.setValue(value);
      return oldValue;
    } else {
      putUnchecked(key, value);
      return null;
    }
  }

  public SWIGTYPE_p_llvm__SparseBitVectorT_t remove(java.lang.Object key) {
    if (!(key instanceof Long)) {
      return null;
    }

    Iterator itr = find((Long) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_llvm__SparseBitVectorT_t oldValue = itr.getValue();
      removeUnchecked(itr);
      return oldValue;
    } else {
      return null;
    }
  }

  public java.util.Set<Entry<Long, SWIGTYPE_p_llvm__SparseBitVectorT_t>> entrySet() {
    java.util.Set<Entry<Long, SWIGTYPE_p_llvm__SparseBitVectorT_t>> setToReturn =
        new java.util.HashSet<Entry<Long, SWIGTYPE_p_llvm__SparseBitVectorT_t>>();

    Iterator itr = begin();
    final Iterator end = end();
    while (itr.isNot(end)) {
      setToReturn.add(new Entry<Long, SWIGTYPE_p_llvm__SparseBitVectorT_t>() {
        private Iterator iterator;

        private Entry<Long, SWIGTYPE_p_llvm__SparseBitVectorT_t> init(Iterator iterator) {
          this.iterator = iterator;
          return this;
        }

        public Long getKey() {
          return iterator.getKey();
        }

        public SWIGTYPE_p_llvm__SparseBitVectorT_t getValue() {
          return iterator.getValue();
        }

        public SWIGTYPE_p_llvm__SparseBitVectorT_t setValue(SWIGTYPE_p_llvm__SparseBitVectorT_t newValue) {
          SWIGTYPE_p_llvm__SparseBitVectorT_t oldValue = iterator.getValue();
          iterator.setValue(newValue);
          return oldValue;
        }
      }.init(itr));
      itr = itr.getNextUnchecked();
    }

    return setToReturn;
  }

  public PtrToBVPtsMap() {
    this(SUPAJNI.new_PtrToBVPtsMap__SWIG_0(), true);
  }

  public PtrToBVPtsMap(PtrToBVPtsMap other) {
    this(SUPAJNI.new_PtrToBVPtsMap__SWIG_1(PtrToBVPtsMap.getCPtr(other), other), true);
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
          SUPAJNI.delete_PtrToBVPtsMap_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator getNextUnchecked() {
      return new SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator(SUPAJNI.PtrToBVPtsMap_Iterator_getNextUnchecked(swigCPtr), true);
    }
  
    private boolean isNot(SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator other) {
      return SUPAJNI.PtrToBVPtsMap_Iterator_isNot(swigCPtr, SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator.getCPtr(other));
    }
  
    private long getKey() {
      return SUPAJNI.PtrToBVPtsMap_Iterator_getKey(swigCPtr);
    }
  
    private SWIGTYPE_p_llvm__SparseBitVectorT_t getValue() {
      return new SWIGTYPE_p_llvm__SparseBitVectorT_t(SUPAJNI.PtrToBVPtsMap_Iterator_getValue(swigCPtr), true);
    }
  
    private void setValue(SWIGTYPE_p_llvm__SparseBitVectorT_t newValue) {
      SUPAJNI.PtrToBVPtsMap_Iterator_setValue(swigCPtr, SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(newValue));
    }
  
  }

  public boolean isEmpty() {
    return SUPAJNI.PtrToBVPtsMap_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.PtrToBVPtsMap_clear(swigCPtr, this);
  }

  private SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator find(long key) {
    return new SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator(SUPAJNI.PtrToBVPtsMap_find(swigCPtr, this, key), true);
  }

  private SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator begin() {
    return new SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator(SUPAJNI.PtrToBVPtsMap_begin(swigCPtr, this), true);
  }

  private SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator end() {
    return new SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator(SUPAJNI.PtrToBVPtsMap_end(swigCPtr, this), true);
  }

  private int sizeImpl() {
    return SUPAJNI.PtrToBVPtsMap_sizeImpl(swigCPtr, this);
  }

  private boolean containsImpl(long key) {
    return SUPAJNI.PtrToBVPtsMap_containsImpl(swigCPtr, this, key);
  }

  private void putUnchecked(long key, SWIGTYPE_p_llvm__SparseBitVectorT_t value) {
    SUPAJNI.PtrToBVPtsMap_putUnchecked(swigCPtr, this, key, SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(value));
  }

  private void removeUnchecked(SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator itr) {
    SUPAJNI.PtrToBVPtsMap_removeUnchecked(swigCPtr, this, SWIGTYPE_p_std__mapT_unsigned_int_llvm__SparseBitVectorT_t_std__lessT_unsigned_int_t_t__iterator.getCPtr(itr));
  }

}
