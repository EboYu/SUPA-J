/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class GlobalDefToRepMapTy extends java.util.AbstractMap<SWIGTYPE_p_llvm__GlobalVariable, SWIGTYPE_p_llvm__GlobalVariable> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected GlobalDefToRepMapTy(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(GlobalDefToRepMapTy obj) {
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
        SUPAJNI.delete_GlobalDefToRepMapTy(swigCPtr);
      }
      swigCPtr = 0;
    }
  }


  public int size() {
    return sizeImpl();
  }

  public boolean containsKey(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_llvm__GlobalVariable)) {
      return false;
    }

    return containsImpl((SWIGTYPE_p_llvm__GlobalVariable)key);
  }

  public SWIGTYPE_p_llvm__GlobalVariable get(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_llvm__GlobalVariable)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_llvm__GlobalVariable) key);
    if (itr.isNot(end())) {
      return itr.getValue();
    }

    return null;
  }

  public SWIGTYPE_p_llvm__GlobalVariable put(SWIGTYPE_p_llvm__GlobalVariable key, SWIGTYPE_p_llvm__GlobalVariable value) {
    Iterator itr = find((SWIGTYPE_p_llvm__GlobalVariable) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_llvm__GlobalVariable oldValue = itr.getValue();
      itr.setValue(value);
      return oldValue;
    } else {
      putUnchecked(key, value);
      return null;
    }
  }

  public SWIGTYPE_p_llvm__GlobalVariable remove(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_llvm__GlobalVariable)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_llvm__GlobalVariable) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_llvm__GlobalVariable oldValue = itr.getValue();
      removeUnchecked(itr);
      return oldValue;
    } else {
      return null;
    }
  }

  public java.util.Set<Entry<SWIGTYPE_p_llvm__GlobalVariable, SWIGTYPE_p_llvm__GlobalVariable>> entrySet() {
    java.util.Set<Entry<SWIGTYPE_p_llvm__GlobalVariable, SWIGTYPE_p_llvm__GlobalVariable>> setToReturn =
        new java.util.HashSet<Entry<SWIGTYPE_p_llvm__GlobalVariable, SWIGTYPE_p_llvm__GlobalVariable>>();

    Iterator itr = begin();
    final Iterator end = end();
    while (itr.isNot(end)) {
      setToReturn.add(new Entry<SWIGTYPE_p_llvm__GlobalVariable, SWIGTYPE_p_llvm__GlobalVariable>() {
        private Iterator iterator;

        private Entry<SWIGTYPE_p_llvm__GlobalVariable, SWIGTYPE_p_llvm__GlobalVariable> init(Iterator iterator) {
          this.iterator = iterator;
          return this;
        }

        public SWIGTYPE_p_llvm__GlobalVariable getKey() {
          return iterator.getKey();
        }

        public SWIGTYPE_p_llvm__GlobalVariable getValue() {
          return iterator.getValue();
        }

        public SWIGTYPE_p_llvm__GlobalVariable setValue(SWIGTYPE_p_llvm__GlobalVariable newValue) {
          SWIGTYPE_p_llvm__GlobalVariable oldValue = iterator.getValue();
          iterator.setValue(newValue);
          return oldValue;
        }
      }.init(itr));
      itr = itr.getNextUnchecked();
    }

    return setToReturn;
  }

  public GlobalDefToRepMapTy() {
    this(SUPAJNI.new_GlobalDefToRepMapTy__SWIG_0(), true);
  }

  public GlobalDefToRepMapTy(GlobalDefToRepMapTy other) {
    this(SUPAJNI.new_GlobalDefToRepMapTy__SWIG_1(GlobalDefToRepMapTy.getCPtr(other), other), true);
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
          SUPAJNI.delete_GlobalDefToRepMapTy_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator getNextUnchecked() {
      return new SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator(SUPAJNI.GlobalDefToRepMapTy_Iterator_getNextUnchecked(swigCPtr), true);
    }
  
    private boolean isNot(SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator other) {
      return SUPAJNI.GlobalDefToRepMapTy_Iterator_isNot(swigCPtr, SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator.getCPtr(other));
    }
  
    private SWIGTYPE_p_llvm__GlobalVariable getKey() {
      long cPtr = SUPAJNI.GlobalDefToRepMapTy_Iterator_getKey(swigCPtr);
      return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__GlobalVariable(cPtr, false);
    }
  
    private SWIGTYPE_p_llvm__GlobalVariable getValue() {
      long cPtr = SUPAJNI.GlobalDefToRepMapTy_Iterator_getValue(swigCPtr);
      return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__GlobalVariable(cPtr, false);
    }
  
    private void setValue(SWIGTYPE_p_llvm__GlobalVariable newValue) {
      SUPAJNI.GlobalDefToRepMapTy_Iterator_setValue(swigCPtr, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(newValue));
    }
  
  }

  public boolean isEmpty() {
    return SUPAJNI.GlobalDefToRepMapTy_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.GlobalDefToRepMapTy_clear(swigCPtr, this);
  }

  private SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator find(SWIGTYPE_p_llvm__GlobalVariable key) {
    return new SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator(SUPAJNI.GlobalDefToRepMapTy_find(swigCPtr, this, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(key)), true);
  }

  private SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator begin() {
    return new SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator(SUPAJNI.GlobalDefToRepMapTy_begin(swigCPtr, this), true);
  }

  private SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator end() {
    return new SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator(SUPAJNI.GlobalDefToRepMapTy_end(swigCPtr, this), true);
  }

  private int sizeImpl() {
    return SUPAJNI.GlobalDefToRepMapTy_sizeImpl(swigCPtr, this);
  }

  private boolean containsImpl(SWIGTYPE_p_llvm__GlobalVariable key) {
    return SUPAJNI.GlobalDefToRepMapTy_containsImpl(swigCPtr, this, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(key));
  }

  private void putUnchecked(SWIGTYPE_p_llvm__GlobalVariable key, SWIGTYPE_p_llvm__GlobalVariable value) {
    SUPAJNI.GlobalDefToRepMapTy_putUnchecked(swigCPtr, this, SWIGTYPE_p_llvm__GlobalVariable.getCPtr(key), SWIGTYPE_p_llvm__GlobalVariable.getCPtr(value));
  }

  private void removeUnchecked(SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator itr) {
    SUPAJNI.GlobalDefToRepMapTy_removeUnchecked(swigCPtr, this, SWIGTYPE_p_std__mapT_llvm__GlobalVariable_const_p_llvm__GlobalVariable_p_std__lessT_llvm__GlobalVariable_const_p_t_t__iterator.getCPtr(itr));
  }

}