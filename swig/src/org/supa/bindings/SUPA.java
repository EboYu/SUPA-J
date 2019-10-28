/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class SUPA implements SUPAConstants {
  public static long imaxabs(long n) {
    return SUPAJNI.imaxabs(n);
  }

  public static imaxdiv_t imaxdiv(long numer, long denom) {
    return new imaxdiv_t(SUPAJNI.imaxdiv(numer, denom), true);
  }

  public static long strtoimax(String nptr, SWIGTYPE_p_p_char endptr, int base) {
    return SUPAJNI.strtoimax(nptr, SWIGTYPE_p_p_char.getCPtr(endptr), base);
  }

  public static java.math.BigInteger strtoumax(String nptr, SWIGTYPE_p_p_char endptr, int base) {
    return SUPAJNI.strtoumax(nptr, SWIGTYPE_p_p_char.getCPtr(endptr), base);
  }

  public static byte[] cdata(SWIGTYPE_p_void ptr, int nelements) {
    return SUPAJNI.cdata(SWIGTYPE_p_void.getCPtr(ptr), nelements);
  }

  public static void memmove(SWIGTYPE_p_void data, byte[] indata) {
    SUPAJNI.memmove(SWIGTYPE_p_void.getCPtr(data), indata);
  }

  public static double cos(double x) {
    return SUPAJNI.cos(x);
  }

  public static double sin(double x) {
    return SUPAJNI.sin(x);
  }

  public static double tan(double x) {
    return SUPAJNI.tan(x);
  }

  public static double acos(double x) {
    return SUPAJNI.acos(x);
  }

  public static double asin(double x) {
    return SUPAJNI.asin(x);
  }

  public static double atan(double x) {
    return SUPAJNI.atan(x);
  }

  public static double atan2(double y, double x) {
    return SUPAJNI.atan2(y, x);
  }

  public static double cosh(double x) {
    return SUPAJNI.cosh(x);
  }

  public static double sinh(double x) {
    return SUPAJNI.sinh(x);
  }

  public static double tanh(double x) {
    return SUPAJNI.tanh(x);
  }

  public static double exp(double x) {
    return SUPAJNI.exp(x);
  }

  public static double log(double x) {
    return SUPAJNI.log(x);
  }

  public static double log10(double x) {
    return SUPAJNI.log10(x);
  }

  public static double pow(double x, double y) {
    return SUPAJNI.pow(x, y);
  }

  public static double sqrt(double x) {
    return SUPAJNI.sqrt(x);
  }

  public static double fabs(double x) {
    return SUPAJNI.fabs(x);
  }

  public static double ceil(double x) {
    return SUPAJNI.ceil(x);
  }

  public static double floor(double x) {
    return SUPAJNI.floor(x);
  }

  public static double fmod(double x, double y) {
    return SUPAJNI.fmod(x, y);
  }

  public static SWIGTYPE_p_llvm__raw_ostream outs() {
    return new SWIGTYPE_p_llvm__raw_ostream(SUPAJNI.outs(), false);
  }

  public static SWIGTYPE_p_llvm__raw_ostream errs() {
    return new SWIGTYPE_p_llvm__raw_ostream(SUPAJNI.errs(), false);
  }

  public static boolean isIntrinsicDbgFun(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isIntrinsicDbgFun(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isInstrinsicDbgInst(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isInstrinsicDbgInst(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isAnAllocationWraper(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isAnAllocationWraper(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isCallSite(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isCallSite(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isNonInstricCallSite(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isNonInstricCallSite(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isReturn(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isReturn(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static SWIGTYPE_p_llvm__Function getLLVMFunction(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.getLLVMFunction(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__CallSite getLLVMCallSite(SWIGTYPE_p_llvm__Instruction inst) {
    return new SWIGTYPE_p_llvm__CallSite(SUPAJNI.getLLVMCallSite(SWIGTYPE_p_llvm__Instruction.getCPtr(inst)), true);
  }

  public static SWIGTYPE_p_llvm__Function getDefFunForMultipleModule(SWIGTYPE_p_llvm__Function fun) {
    long cPtr = SUPAJNI.getDefFunForMultipleModule(SWIGTYPE_p_llvm__Function.getCPtr(fun));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Function getCallee(SWIGTYPE_p_llvm__CallSite cs) {
    long cPtr = SUPAJNI.getCallee__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Function getCallee(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getCallee__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
  }

  public static boolean isExtCall(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isExtCall__SWIG_0(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isExtCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isExtCall__SWIG_1(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isExtCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isExtCall__SWIG_2(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isHeapAllocExtFunViaRet(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isHeapAllocExtFunViaRet(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isHeapAllocExtFunViaArg(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isHeapAllocExtFunViaArg(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isHeapAllocExtCallViaRet(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isHeapAllocExtCallViaRet__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isHeapAllocExtCallViaRet(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isHeapAllocExtCallViaRet__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isHeapAllocExtCallViaArg(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isHeapAllocExtCallViaArg__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isHeapAllocExtCallViaArg(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isHeapAllocExtCallViaArg__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isHeapAllocExtCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isHeapAllocExtCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isHeapAllocExtCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isHeapAllocExtCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static int getHeapAllocHoldingArgPosition(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.getHeapAllocHoldingArgPosition__SWIG_0(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static int getHeapAllocHoldingArgPosition(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.getHeapAllocHoldingArgPosition__SWIG_1(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static int getHeapAllocHoldingArgPosition(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.getHeapAllocHoldingArgPosition__SWIG_2(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isReallocExtFun(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isReallocExtFun(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isReallocExtCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isReallocExtCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isReallocExtCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isReallocExtCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isDeallocExtFun(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isDeallocExtFun(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isDeallocExtCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isDeallocExtCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isDeallocExtCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isDeallocExtCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isStaticExtFun(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isStaticExtFun(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isStaticExtCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isStaticExtCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isStaticExtCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isStaticExtCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isHeapAllocOrStaticExtCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isHeapAllocOrStaticExtCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isHeapAllocOrStaticExtCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isHeapAllocOrStaticExtCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static SWIGTYPE_p_ExtAPI__extf_t extCallTy(SWIGTYPE_p_llvm__Function fun) {
    return new SWIGTYPE_p_ExtAPI__extf_t(SUPAJNI.extCallTy(SWIGTYPE_p_llvm__Function.getCPtr(fun)), true);
  }

  public static SWIGTYPE_p_llvm__PointerType getRefTypeOfHeapAllocOrStatic(SWIGTYPE_p_llvm__CallSite cs) {
    long cPtr = SUPAJNI.getRefTypeOfHeapAllocOrStatic__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__PointerType(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__PointerType getRefTypeOfHeapAllocOrStatic(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getRefTypeOfHeapAllocOrStatic__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__PointerType(cPtr, false);
  }

  public static boolean isThreadForkCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isThreadForkCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isThreadForkCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isThreadForkCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isHareParForCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isHareParForCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isHareParForCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isHareParForCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isThreadJoinCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isThreadJoinCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isThreadJoinCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isThreadJoinCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isThreadExitCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isThreadExitCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isThreadExitCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isThreadExitCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isLockAquireCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isLockAquireCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isLockAquireCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isLockAquireCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isLockReleaseCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isLockReleaseCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isLockReleaseCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isLockReleaseCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isBarrierWaitCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isBarrierWaitCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isBarrierWaitCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isBarrierWaitCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static SWIGTYPE_p_llvm__Value getForkedFun(SWIGTYPE_p_llvm__CallSite cs) {
    long cPtr = SUPAJNI.getForkedFun__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getForkedFun(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getForkedFun__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getActualParmAtForkSite(SWIGTYPE_p_llvm__CallSite cs) {
    long cPtr = SUPAJNI.getActualParmAtForkSite__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getActualParmAtForkSite(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getActualParmAtForkSite__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getTaskFuncAtHareParForSite(SWIGTYPE_p_llvm__CallSite cs) {
    long cPtr = SUPAJNI.getTaskFuncAtHareParForSite__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getTaskFuncAtHareParForSite(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getTaskFuncAtHareParForSite__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getTaskDataAtHareParForSite(SWIGTYPE_p_llvm__CallSite cs) {
    long cPtr = SUPAJNI.getTaskDataAtHareParForSite__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value getTaskDataAtHareParForSite(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getTaskDataAtHareParForSite__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static boolean isObject(SWIGTYPE_p_llvm__Value ref) {
    return SUPAJNI.isObject(SWIGTYPE_p_llvm__Value.getCPtr(ref));
  }

  public static boolean isDeadFunction(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isDeadFunction(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean ArgInDeadFunction(SWIGTYPE_p_llvm__Value val) {
    return SUPAJNI.ArgInDeadFunction(SWIGTYPE_p_llvm__Value.getCPtr(val));
  }

  public static boolean isProgEntryFunction(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isProgEntryFunction(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static SWIGTYPE_p_llvm__Function getProgEntryFunction(SVFModule svfModule) {
    long cPtr = SUPAJNI.getProgEntryFunction(SVFModule.getCPtr(svfModule), svfModule);
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Function(cPtr, false);
  }

  public static boolean ArgInProgEntryFunction(SWIGTYPE_p_llvm__Value val) {
    return SUPAJNI.ArgInProgEntryFunction(SWIGTYPE_p_llvm__Value.getCPtr(val));
  }

  public static boolean isPtrInDeadFunction(SWIGTYPE_p_llvm__Value value) {
    return SUPAJNI.isPtrInDeadFunction(SWIGTYPE_p_llvm__Value.getCPtr(value));
  }

  public static boolean isProgExitFunction(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isProgExitFunction(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean isProgExitCall(SWIGTYPE_p_llvm__CallSite cs) {
    return SUPAJNI.isProgExitCall__SWIG_0(SWIGTYPE_p_llvm__CallSite.getCPtr(cs));
  }

  public static boolean isProgExitCall(SWIGTYPE_p_llvm__Instruction inst) {
    return SUPAJNI.isProgExitCall__SWIG_1(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
  }

  public static boolean isNoCallerFunction(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.isNoCallerFunction(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static boolean ArgInNoCallerFunction(SWIGTYPE_p_llvm__Value val) {
    return SUPAJNI.ArgInNoCallerFunction(SWIGTYPE_p_llvm__Value.getCPtr(val));
  }

  public static boolean functionDoesNotRet(SWIGTYPE_p_llvm__Function fun) {
    return SUPAJNI.functionDoesNotRet(SWIGTYPE_p_llvm__Function.getCPtr(fun));
  }

  public static void getFunReachableBBs(SWIGTYPE_p_llvm__Function fun, SWIGTYPE_p_llvm__DominatorTree dt, BasicBlockVector bbs) {
    SUPAJNI.getFunReachableBBs(SWIGTYPE_p_llvm__Function.getCPtr(fun), SWIGTYPE_p_llvm__DominatorTree.getCPtr(dt), BasicBlockVector.getCPtr(bbs), bbs);
  }

  public static SWIGTYPE_p_llvm__BasicBlock getFunExitBB(SWIGTYPE_p_llvm__Function fun) {
    long cPtr = SUPAJNI.getFunExitBB(SWIGTYPE_p_llvm__Function.getCPtr(fun));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__BasicBlock(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value stripConstantCasts(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.stripConstantCasts(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Value stripAllCasts(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.stripAllCasts(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Value(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__Type getTypeOfHeapAlloc(SWIGTYPE_p_llvm__Instruction inst) {
    long cPtr = SUPAJNI.getTypeOfHeapAlloc(SWIGTYPE_p_llvm__Instruction.getCPtr(inst));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__Type(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isGepConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isGepConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isInt2PtrConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isInt2PtrConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isPtr2IntConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isPtr2IntConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isCastConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isCastConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isSelectConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isSelectConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isTruncConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isTruncConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static SWIGTYPE_p_llvm__ConstantExpr isCmpConstantExpr(SWIGTYPE_p_llvm__Value val) {
    long cPtr = SUPAJNI.isCmpConstantExpr(SWIGTYPE_p_llvm__Value.getCPtr(val));
    return (cPtr == 0) ? null : new SWIGTYPE_p_llvm__ConstantExpr(cPtr, false);
  }

  public static void getNextInsts(SWIGTYPE_p_llvm__Instruction curInst, InstructionVector instList) {
    SUPAJNI.getNextInsts(SWIGTYPE_p_llvm__Instruction.getCPtr(curInst), InstructionVector.getCPtr(instList), instList);
  }

  public static void getPrevInsts(SWIGTYPE_p_llvm__Instruction curInst, InstructionVector instList) {
    SUPAJNI.getPrevInsts(SWIGTYPE_p_llvm__Instruction.getCPtr(curInst), InstructionVector.getCPtr(instList), instList);
  }

  public static long getBBSuccessorPos(SWIGTYPE_p_llvm__BasicBlock BB, SWIGTYPE_p_llvm__BasicBlock Succ) {
    return SUPAJNI.getBBSuccessorPos(SWIGTYPE_p_llvm__BasicBlock.getCPtr(BB), SWIGTYPE_p_llvm__BasicBlock.getCPtr(Succ));
  }

  public static long getBBSuccessorNum(SWIGTYPE_p_llvm__BasicBlock BB) {
    return SUPAJNI.getBBSuccessorNum(SWIGTYPE_p_llvm__BasicBlock.getCPtr(BB));
  }

  public static long getBBPredecessorPos(SWIGTYPE_p_llvm__BasicBlock BB, SWIGTYPE_p_llvm__BasicBlock Pred) {
    return SUPAJNI.getBBPredecessorPos(SWIGTYPE_p_llvm__BasicBlock.getCPtr(BB), SWIGTYPE_p_llvm__BasicBlock.getCPtr(Pred));
  }

  public static long getBBPredecessorNum(SWIGTYPE_p_llvm__BasicBlock BB) {
    return SUPAJNI.getBBPredecessorNum(SWIGTYPE_p_llvm__BasicBlock.getCPtr(BB));
  }

  public static String getSourceLoc(SWIGTYPE_p_llvm__Value val) {
    return SUPAJNI.getSourceLoc(SWIGTYPE_p_llvm__Value.getCPtr(val));
  }

  public static String getSourceLocOfFunction(SWIGTYPE_p_llvm__Function F) {
    return SUPAJNI.getSourceLocOfFunction(SWIGTYPE_p_llvm__Function.getCPtr(F));
  }

  public static void dumpSet(SWIGTYPE_p_llvm__SparseBitVectorT_t To, SWIGTYPE_p_llvm__raw_ostream O) {
    SUPAJNI.dumpSet__SWIG_0(SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(To), SWIGTYPE_p_llvm__raw_ostream.getCPtr(O));
  }

  public static void dumpSet(SWIGTYPE_p_llvm__SparseBitVectorT_t To) {
    SUPAJNI.dumpSet__SWIG_1(SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(To));
  }

  public static void dumpPointsToSet(long node, SWIGTYPE_p_llvm__SparseBitVectorT_t To) {
    SUPAJNI.dumpPointsToSet(node, SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(To));
  }

  public static void dumpAliasSet(long node, SWIGTYPE_p_llvm__SparseBitVectorT_t To) {
    SUPAJNI.dumpAliasSet(node, SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(To));
  }

  public static String sucMsg(String msg) {
    return SUPAJNI.sucMsg(msg);
  }

  public static void wrnMsg(String msg) {
    SUPAJNI.wrnMsg(msg);
  }

  public static String errMsg(String msg) {
    return SUPAJNI.errMsg(msg);
  }

  public static String bugMsg1(String msg) {
    return SUPAJNI.bugMsg1(msg);
  }

  public static String bugMsg2(String msg) {
    return SUPAJNI.bugMsg2(msg);
  }

  public static String bugMsg3(String msg) {
    return SUPAJNI.bugMsg3(msg);
  }

  public static String pasMsg(String msg) {
    return SUPAJNI.pasMsg(msg);
  }

  public static void reportMemoryUsageKB(String infor, SWIGTYPE_p_llvm__raw_ostream O) {
    SUPAJNI.reportMemoryUsageKB__SWIG_0(infor, SWIGTYPE_p_llvm__raw_ostream.getCPtr(O));
  }

  public static void reportMemoryUsageKB(String infor) {
    SUPAJNI.reportMemoryUsageKB__SWIG_1(infor);
  }

  public static boolean getMemoryUsageKB(SWIGTYPE_p_unsigned_int vmrss_kb, SWIGTYPE_p_unsigned_int vmsize_kb) {
    return SUPAJNI.getMemoryUsageKB(SWIGTYPE_p_unsigned_int.getCPtr(vmrss_kb), SWIGTYPE_p_unsigned_int.getCPtr(vmsize_kb));
  }

  public static void increaseStackSize() {
    SUPAJNI.increaseStackSize();
  }

  public static boolean isIRFile(String filename) {
    return SUPAJNI.isIRFile(filename);
  }

  public static void processArguments(int argc, SWIGTYPE_p_p_char argv, SWIGTYPE_p_int arg_num, SWIGTYPE_p_p_char arg_value, StringVector moduleNameVec) {
    SUPAJNI.processArguments(argc, SWIGTYPE_p_p_char.getCPtr(argv), SWIGTYPE_p_int.getCPtr(arg_num), SWIGTYPE_p_p_char.getCPtr(arg_value), StringVector.getCPtr(moduleNameVec), moduleNameVec);
  }

  public static boolean cmpPts(SWIGTYPE_p_llvm__SparseBitVectorT_t lpts, SWIGTYPE_p_llvm__SparseBitVectorT_t rpts) {
    return SUPAJNI.cmpPts(SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(lpts), SWIGTYPE_p_llvm__SparseBitVectorT_t.getCPtr(rpts));
  }

}