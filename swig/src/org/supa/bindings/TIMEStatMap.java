/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class TIMEStatMap extends java.util.AbstractMap<String, Double> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected TIMEStatMap(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(TIMEStatMap obj) {
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
        SUPAJNI.delete_TIMEStatMap(swigCPtr);
      }
      swigCPtr = 0;
    }
  }


  public int size() {
    return sizeImpl();
  }

  public boolean containsKey(java.lang.Object key) {
    if (!(key instanceof String)) {
      return false;
    }

    return containsImpl((String)key);
  }

  public Double get(java.lang.Object key) {
    if (!(key instanceof String)) {
      return null;
    }

    Iterator itr = find((String) key);
    if (itr.isNot(end())) {
      return itr.getValue();
    }

    return null;
  }

  public Double put(String key, Double value) {
    Iterator itr = find((String) key);
    if (itr.isNot(end())) {
      Double oldValue = itr.getValue();
      itr.setValue(value);
      return oldValue;
    } else {
      putUnchecked(key, value);
      return null;
    }
  }

  public Double remove(java.lang.Object key) {
    if (!(key instanceof String)) {
      return null;
    }

    Iterator itr = find((String) key);
    if (itr.isNot(end())) {
      Double oldValue = itr.getValue();
      removeUnchecked(itr);
      return oldValue;
    } else {
      return null;
    }
  }

  public java.util.Set<Entry<String, Double>> entrySet() {
    java.util.Set<Entry<String, Double>> setToReturn =
        new java.util.HashSet<Entry<String, Double>>();

    Iterator itr = begin();
    final Iterator end = end();
    while (itr.isNot(end)) {
      setToReturn.add(new Entry<String, Double>() {
        private Iterator iterator;

        private Entry<String, Double> init(Iterator iterator) {
          this.iterator = iterator;
          return this;
        }

        public String getKey() {
          return iterator.getKey();
        }

        public Double getValue() {
          return iterator.getValue();
        }

        public Double setValue(Double newValue) {
          Double oldValue = iterator.getValue();
          iterator.setValue(newValue);
          return oldValue;
        }
      }.init(itr));
      itr = itr.getNextUnchecked();
    }

    return setToReturn;
  }

  public TIMEStatMap() {
    this(SUPAJNI.new_TIMEStatMap__SWIG_0(), true);
  }

  public TIMEStatMap(TIMEStatMap other) {
    this(SUPAJNI.new_TIMEStatMap__SWIG_1(TIMEStatMap.getCPtr(other), other), true);
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
          SUPAJNI.delete_TIMEStatMap_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private TIMEStatMap.Iterator getNextUnchecked() {
      return new TIMEStatMap.Iterator(SUPAJNI.TIMEStatMap_Iterator_getNextUnchecked(swigCPtr, this), true);
    }
  
    private boolean isNot(TIMEStatMap.Iterator other) {
      return SUPAJNI.TIMEStatMap_Iterator_isNot(swigCPtr, this, TIMEStatMap.Iterator.getCPtr(other), other);
    }
  
    private String getKey() {
      return SUPAJNI.TIMEStatMap_Iterator_getKey(swigCPtr, this);
    }
  
    private double getValue() {
      return SUPAJNI.TIMEStatMap_Iterator_getValue(swigCPtr, this);
    }
  
    private void setValue(double newValue) {
      SUPAJNI.TIMEStatMap_Iterator_setValue(swigCPtr, this, newValue);
    }
  
  }

  public boolean isEmpty() {
    return SUPAJNI.TIMEStatMap_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.TIMEStatMap_clear(swigCPtr, this);
  }

  private TIMEStatMap.Iterator find(String key) {
    return new TIMEStatMap.Iterator(SUPAJNI.TIMEStatMap_find(swigCPtr, this, key), true);
  }

  private TIMEStatMap.Iterator begin() {
    return new TIMEStatMap.Iterator(SUPAJNI.TIMEStatMap_begin(swigCPtr, this), true);
  }

  private TIMEStatMap.Iterator end() {
    return new TIMEStatMap.Iterator(SUPAJNI.TIMEStatMap_end(swigCPtr, this), true);
  }

  private int sizeImpl() {
    return SUPAJNI.TIMEStatMap_sizeImpl(swigCPtr, this);
  }

  private boolean containsImpl(String key) {
    return SUPAJNI.TIMEStatMap_containsImpl(swigCPtr, this, key);
  }

  private void putUnchecked(String key, double value) {
    SUPAJNI.TIMEStatMap_putUnchecked(swigCPtr, this, key, value);
  }

  private void removeUnchecked(TIMEStatMap.Iterator itr) {
    SUPAJNI.TIMEStatMap_removeUnchecked(swigCPtr, this, TIMEStatMap.Iterator.getCPtr(itr), itr);
  }

}
