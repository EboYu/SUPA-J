/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class DPMToDPMMap extends java.util.AbstractMap<SWIGTYPE_p_DPItem, SWIGTYPE_p_DPItem> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected DPMToDPMMap(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(DPMToDPMMap obj) {
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
        SUPAJNI.delete_DPMToDPMMap(swigCPtr);
      }
      swigCPtr = 0;
    }
  }


  public int size() {
    return sizeImpl();
  }

  public boolean containsKey(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_DPItem)) {
      return false;
    }

    return containsImpl((SWIGTYPE_p_DPItem)key);
  }

  public SWIGTYPE_p_DPItem get(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_DPItem)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_DPItem) key);
    if (itr.isNot(end())) {
      return itr.getValue();
    }

    return null;
  }

  public SWIGTYPE_p_DPItem put(SWIGTYPE_p_DPItem key, SWIGTYPE_p_DPItem value) {
    Iterator itr = find((SWIGTYPE_p_DPItem) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_DPItem oldValue = itr.getValue();
      itr.setValue(value);
      return oldValue;
    } else {
      putUnchecked(key, value);
      return null;
    }
  }

  public SWIGTYPE_p_DPItem remove(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_DPItem)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_DPItem) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_DPItem oldValue = itr.getValue();
      removeUnchecked(itr);
      return oldValue;
    } else {
      return null;
    }
  }

  public java.util.Set<Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_DPItem>> entrySet() {
    java.util.Set<Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_DPItem>> setToReturn =
        new java.util.HashSet<Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_DPItem>>();

    Iterator itr = begin();
    final Iterator end = end();
    while (itr.isNot(end)) {
      setToReturn.add(new Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_DPItem>() {
        private Iterator iterator;

        private Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_DPItem> init(Iterator iterator) {
          this.iterator = iterator;
          return this;
        }

        public SWIGTYPE_p_DPItem getKey() {
          return iterator.getKey();
        }

        public SWIGTYPE_p_DPItem getValue() {
          return iterator.getValue();
        }

        public SWIGTYPE_p_DPItem setValue(SWIGTYPE_p_DPItem newValue) {
          SWIGTYPE_p_DPItem oldValue = iterator.getValue();
          iterator.setValue(newValue);
          return oldValue;
        }
      }.init(itr));
      itr = itr.getNextUnchecked();
    }

    return setToReturn;
  }

  public DPMToDPMMap() {
    this(SUPAJNI.new_DPMToDPMMap__SWIG_0(), true);
  }

  public DPMToDPMMap(DPMToDPMMap other) {
    this(SUPAJNI.new_DPMToDPMMap__SWIG_1(DPMToDPMMap.getCPtr(other), other), true);
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
          SUPAJNI.delete_DPMToDPMMap_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private DPMToDPMMap.Iterator getNextUnchecked() {
      return new DPMToDPMMap.Iterator(SUPAJNI.DPMToDPMMap_Iterator_getNextUnchecked(swigCPtr, this), true);
    }
  
    private boolean isNot(DPMToDPMMap.Iterator other) {
      return SUPAJNI.DPMToDPMMap_Iterator_isNot(swigCPtr, this, DPMToDPMMap.Iterator.getCPtr(other), other);
    }
  
    private SWIGTYPE_p_DPItem getKey() {
      return new SWIGTYPE_p_DPItem(SUPAJNI.DPMToDPMMap_Iterator_getKey(swigCPtr, this), true);
    }
  
    private SWIGTYPE_p_DPItem getValue() {
      return new SWIGTYPE_p_DPItem(SUPAJNI.DPMToDPMMap_Iterator_getValue(swigCPtr, this), true);
    }
  
    private void setValue(SWIGTYPE_p_DPItem newValue) {
      SUPAJNI.DPMToDPMMap_Iterator_setValue(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(newValue));
    }
  
  }

  public boolean isEmpty() {
    return SUPAJNI.DPMToDPMMap_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.DPMToDPMMap_clear(swigCPtr, this);
  }

  private DPMToDPMMap.Iterator find(SWIGTYPE_p_DPItem key) {
    return new DPMToDPMMap.Iterator(SUPAJNI.DPMToDPMMap_find(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(key)), true);
  }

  private DPMToDPMMap.Iterator begin() {
    return new DPMToDPMMap.Iterator(SUPAJNI.DPMToDPMMap_begin(swigCPtr, this), true);
  }

  private DPMToDPMMap.Iterator end() {
    return new DPMToDPMMap.Iterator(SUPAJNI.DPMToDPMMap_end(swigCPtr, this), true);
  }

  private int sizeImpl() {
    return SUPAJNI.DPMToDPMMap_sizeImpl(swigCPtr, this);
  }

  private boolean containsImpl(SWIGTYPE_p_DPItem key) {
    return SUPAJNI.DPMToDPMMap_containsImpl(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(key));
  }

  private void putUnchecked(SWIGTYPE_p_DPItem key, SWIGTYPE_p_DPItem value) {
    SUPAJNI.DPMToDPMMap_putUnchecked(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(key), SWIGTYPE_p_DPItem.getCPtr(value));
  }

  private void removeUnchecked(DPMToDPMMap.Iterator itr) {
    SUPAJNI.DPMToDPMMap_removeUnchecked(swigCPtr, this, DPMToDPMMap.Iterator.getCPtr(itr), itr);
  }

}