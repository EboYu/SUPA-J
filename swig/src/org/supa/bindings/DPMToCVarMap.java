/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class DPMToCVarMap extends java.util.AbstractMap<SWIGTYPE_p_DPItem, SWIGTYPE_p_CVar> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected DPMToCVarMap(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(DPMToCVarMap obj) {
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
        SUPAJNI.delete_DPMToCVarMap(swigCPtr);
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

  public SWIGTYPE_p_CVar get(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_DPItem)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_DPItem) key);
    if (itr.isNot(end())) {
      return itr.getValue();
    }

    return null;
  }

  public SWIGTYPE_p_CVar put(SWIGTYPE_p_DPItem key, SWIGTYPE_p_CVar value) {
    Iterator itr = find((SWIGTYPE_p_DPItem) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_CVar oldValue = itr.getValue();
      itr.setValue(value);
      return oldValue;
    } else {
      putUnchecked(key, value);
      return null;
    }
  }

  public SWIGTYPE_p_CVar remove(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_DPItem)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_DPItem) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_CVar oldValue = itr.getValue();
      removeUnchecked(itr);
      return oldValue;
    } else {
      return null;
    }
  }

  public java.util.Set<Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_CVar>> entrySet() {
    java.util.Set<Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_CVar>> setToReturn =
        new java.util.HashSet<Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_CVar>>();

    Iterator itr = begin();
    final Iterator end = end();
    while (itr.isNot(end)) {
      setToReturn.add(new Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_CVar>() {
        private Iterator iterator;

        private Entry<SWIGTYPE_p_DPItem, SWIGTYPE_p_CVar> init(Iterator iterator) {
          this.iterator = iterator;
          return this;
        }

        public SWIGTYPE_p_DPItem getKey() {
          return iterator.getKey();
        }

        public SWIGTYPE_p_CVar getValue() {
          return iterator.getValue();
        }

        public SWIGTYPE_p_CVar setValue(SWIGTYPE_p_CVar newValue) {
          SWIGTYPE_p_CVar oldValue = iterator.getValue();
          iterator.setValue(newValue);
          return oldValue;
        }
      }.init(itr));
      itr = itr.getNextUnchecked();
    }

    return setToReturn;
  }

  public DPMToCVarMap() {
    this(SUPAJNI.new_DPMToCVarMap__SWIG_0(), true);
  }

  public DPMToCVarMap(DPMToCVarMap other) {
    this(SUPAJNI.new_DPMToCVarMap__SWIG_1(DPMToCVarMap.getCPtr(other), other), true);
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
          SUPAJNI.delete_DPMToCVarMap_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private DPMToCVarMap.Iterator getNextUnchecked() {
      return new DPMToCVarMap.Iterator(SUPAJNI.DPMToCVarMap_Iterator_getNextUnchecked(swigCPtr, this), true);
    }
  
    private boolean isNot(DPMToCVarMap.Iterator other) {
      return SUPAJNI.DPMToCVarMap_Iterator_isNot(swigCPtr, this, DPMToCVarMap.Iterator.getCPtr(other), other);
    }
  
    private SWIGTYPE_p_DPItem getKey() {
      return new SWIGTYPE_p_DPItem(SUPAJNI.DPMToCVarMap_Iterator_getKey(swigCPtr, this), true);
    }
  
    private SWIGTYPE_p_CVar getValue() {
      return new SWIGTYPE_p_CVar(SUPAJNI.DPMToCVarMap_Iterator_getValue(swigCPtr, this), true);
    }
  
    private void setValue(SWIGTYPE_p_CVar newValue) {
      SUPAJNI.DPMToCVarMap_Iterator_setValue(swigCPtr, this, SWIGTYPE_p_CVar.getCPtr(newValue));
    }
  
  }

  public boolean isEmpty() {
    return SUPAJNI.DPMToCVarMap_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.DPMToCVarMap_clear(swigCPtr, this);
  }

  private DPMToCVarMap.Iterator find(SWIGTYPE_p_DPItem key) {
    return new DPMToCVarMap.Iterator(SUPAJNI.DPMToCVarMap_find(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(key)), true);
  }

  private DPMToCVarMap.Iterator begin() {
    return new DPMToCVarMap.Iterator(SUPAJNI.DPMToCVarMap_begin(swigCPtr, this), true);
  }

  private DPMToCVarMap.Iterator end() {
    return new DPMToCVarMap.Iterator(SUPAJNI.DPMToCVarMap_end(swigCPtr, this), true);
  }

  private int sizeImpl() {
    return SUPAJNI.DPMToCVarMap_sizeImpl(swigCPtr, this);
  }

  private boolean containsImpl(SWIGTYPE_p_DPItem key) {
    return SUPAJNI.DPMToCVarMap_containsImpl(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(key));
  }

  private void putUnchecked(SWIGTYPE_p_DPItem key, SWIGTYPE_p_CVar value) {
    SUPAJNI.DPMToCVarMap_putUnchecked(swigCPtr, this, SWIGTYPE_p_DPItem.getCPtr(key), SWIGTYPE_p_CVar.getCPtr(value));
  }

  private void removeUnchecked(DPMToCVarMap.Iterator itr) {
    SUPAJNI.DPMToCVarMap_removeUnchecked(swigCPtr, this, DPMToCVarMap.Iterator.getCPtr(itr), itr);
  }

}