/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class FunDefToDeclsMapTy extends java.util.AbstractMap<SWIGTYPE_p_llvm__Function, SWIGTYPE_p_FunctionSetType> {
  private transient long swigCPtr;
  protected transient boolean swigCMemOwn;

  protected FunDefToDeclsMapTy(long cPtr, boolean cMemoryOwn) {
    swigCMemOwn = cMemoryOwn;
    swigCPtr = cPtr;
  }

  protected static long getCPtr(FunDefToDeclsMapTy obj) {
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
        SUPAJNI.delete_FunDefToDeclsMapTy(swigCPtr);
      }
      swigCPtr = 0;
    }
  }


  public int size() {
    return sizeImpl();
  }

  public boolean containsKey(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_llvm__Function)) {
      return false;
    }

    return containsImpl((SWIGTYPE_p_llvm__Function)key);
  }

  public SWIGTYPE_p_FunctionSetType get(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_llvm__Function)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_llvm__Function) key);
    if (itr.isNot(end())) {
      return itr.getValue();
    }

    return null;
  }

  public SWIGTYPE_p_FunctionSetType put(SWIGTYPE_p_llvm__Function key, SWIGTYPE_p_FunctionSetType value) {
    Iterator itr = find((SWIGTYPE_p_llvm__Function) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_FunctionSetType oldValue = itr.getValue();
      itr.setValue(value);
      return oldValue;
    } else {
      putUnchecked(key, value);
      return null;
    }
  }

  public SWIGTYPE_p_FunctionSetType remove(java.lang.Object key) {
    if (!(key instanceof SWIGTYPE_p_llvm__Function)) {
      return null;
    }

    Iterator itr = find((SWIGTYPE_p_llvm__Function) key);
    if (itr.isNot(end())) {
      SWIGTYPE_p_FunctionSetType oldValue = itr.getValue();
      removeUnchecked(itr);
      return oldValue;
    } else {
      return null;
    }
  }

  public java.util.Set<Entry<SWIGTYPE_p_llvm__Function, SWIGTYPE_p_FunctionSetType>> entrySet() {
    java.util.Set<Entry<SWIGTYPE_p_llvm__Function, SWIGTYPE_p_FunctionSetType>> setToReturn =
        new java.util.HashSet<Entry<SWIGTYPE_p_llvm__Function, SWIGTYPE_p_FunctionSetType>>();

    Iterator itr = begin();
    final Iterator end = end();
    while (itr.isNot(end)) {
      setToReturn.add(new Entry<SWIGTYPE_p_llvm__Function, SWIGTYPE_p_FunctionSetType>() {
        private Iterator iterator;

        private Entry<SWIGTYPE_p_llvm__Function, SWIGTYPE_p_FunctionSetType> init(Iterator iterator) {
          this.iterator = iterator;
          return this;
        }

        public SWIGTYPE_p_llvm__Function getKey() {
          return iterator.getKey();
        }

        public SWIGTYPE_p_FunctionSetType getValue() {
          return iterator.getValue();
        }

        public SWIGTYPE_p_FunctionSetType setValue(SWIGTYPE_p_FunctionSetType newValue) {
          SWIGTYPE_p_FunctionSetType oldValue = iterator.getValue();
          iterator.setValue(newValue);
          return oldValue;
        }
      }.init(itr));
      itr = itr.getNextUnchecked();
    }

    return setToReturn;
  }

  public FunDefToDeclsMapTy() {
    this(SUPAJNI.new_FunDefToDeclsMapTy__SWIG_0(), true);
  }

  public FunDefToDeclsMapTy(FunDefToDeclsMapTy other) {
    this(SUPAJNI.new_FunDefToDeclsMapTy__SWIG_1(FunDefToDeclsMapTy.getCPtr(other), other), true);
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
          SUPAJNI.delete_FunDefToDeclsMapTy_Iterator(swigCPtr);
        }
        swigCPtr = 0;
      }
    }
  
    private SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator getNextUnchecked() {
      return new SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator(SUPAJNI.FunDefToDeclsMapTy_Iterator_getNextUnchecked(swigCPtr), true);
    }
  
    private boolean isNot(SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator other) {
      return SUPAJNI.FunDefToDeclsMapTy_Iterator_isNot(swigCPtr, SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator.getCPtr(other));
    }
  
    private SWIGTYPE_p_llvm__Function getKey() {
      long cPtr = SUPAJNI.FunDefToDeclsMapTy_Iterator_getKey(swigCPtr);
      return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
    }
  
    private SWIGTYPE_p_FunctionSetType getValue() {
      return new SWIGTYPE_p_FunctionSetType(SUPAJNI.FunDefToDeclsMapTy_Iterator_getValue(swigCPtr), true);
    }
  
    private void setValue(SWIGTYPE_p_FunctionSetType newValue) {
      SUPAJNI.FunDefToDeclsMapTy_Iterator_setValue(swigCPtr, SWIGTYPE_p_FunctionSetType.getCPtr(newValue));
    }
  
  }

  public boolean isEmpty() {
    return SUPAJNI.FunDefToDeclsMapTy_isEmpty(swigCPtr, this);
  }

  public void clear() {
    SUPAJNI.FunDefToDeclsMapTy_clear(swigCPtr, this);
  }

  private SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator find(SWIGTYPE_p_llvm__Function key) {
    return new SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator(SUPAJNI.FunDefToDeclsMapTy_find(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(key)), true);
  }

  private SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator begin() {
    return new SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator(SUPAJNI.FunDefToDeclsMapTy_begin(swigCPtr, this), true);
  }

  private SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator end() {
    return new SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator(SUPAJNI.FunDefToDeclsMapTy_end(swigCPtr, this), true);
  }

  private int sizeImpl() {
    return SUPAJNI.FunDefToDeclsMapTy_sizeImpl(swigCPtr, this);
  }

  private boolean containsImpl(SWIGTYPE_p_llvm__Function key) {
    return SUPAJNI.FunDefToDeclsMapTy_containsImpl(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(key));
  }

  private void putUnchecked(SWIGTYPE_p_llvm__Function key, SWIGTYPE_p_FunctionSetType value) {
    SUPAJNI.FunDefToDeclsMapTy_putUnchecked(swigCPtr, this, SWIGTYPE_p_llvm__Function.getCPtr(key), SWIGTYPE_p_FunctionSetType.getCPtr(value));
  }

  private void removeUnchecked(SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator itr) {
    SUPAJNI.FunDefToDeclsMapTy_removeUnchecked(swigCPtr, this, SWIGTYPE_p_std__mapT_llvm__Function_const_p_FunctionSetType_std__lessT_llvm__Function_const_p_t_t__iterator.getCPtr(itr));
  }

}