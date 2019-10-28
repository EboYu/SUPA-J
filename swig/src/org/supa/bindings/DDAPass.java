/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 4.0.2
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */

package org.supa.bindings;

public class DDAPass extends ModulePass {
  private transient long swigCPtr;

  protected DDAPass(long cPtr, boolean cMemoryOwn) {
    super(SUPAJNI.DDAPass_SWIGUpcast(cPtr), cMemoryOwn);
    swigCPtr = cPtr;
  }

  protected static long getCPtr(DDAPass obj) {
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
        SUPAJNI.delete_DDAPass(swigCPtr);
      }
      swigCPtr = 0;
    }
    super.delete();
  }

  public static void setID(char value) {
    SUPAJNI.DDAPass_ID_set(value);
  }

  public static char getID() {
    return SUPAJNI.DDAPass_ID_get();
  }

  public DDAPass() {
    this(SUPAJNI.new_DDAPass(), true);
  }

  public void getAnalysisUsage(SWIGTYPE_p_llvm__AnalysisUsage au) {
    SUPAJNI.DDAPass_getAnalysisUsage(swigCPtr, this, SWIGTYPE_p_llvm__AnalysisUsage.getCPtr(au));
  }

  public SWIGTYPE_p_void getAdjustedAnalysisPointer(SWIGTYPE_p_void id) {
    long cPtr = SUPAJNI.DDAPass_getAdjustedAnalysisPointer(swigCPtr, this, SWIGTYPE_p_void.getCPtr(id));
    return (cPtr == 0) ? null : new SWIGTYPE_p_void(cPtr, false);
  }

  public SWIGTYPE_p_llvm__AliasResult alias(SWIGTYPE_p_llvm__MemoryLocation LocA, SWIGTYPE_p_llvm__MemoryLocation LocB) {
    return new SWIGTYPE_p_llvm__AliasResult(SUPAJNI.DDAPass_alias__SWIG_0(swigCPtr, this, SWIGTYPE_p_llvm__MemoryLocation.getCPtr(LocA), SWIGTYPE_p_llvm__MemoryLocation.getCPtr(LocB)), true);
  }

  public SWIGTYPE_p_llvm__AliasResult alias(SWIGTYPE_p_llvm__Value V1, SWIGTYPE_p_llvm__Value V2) {
    return new SWIGTYPE_p_llvm__AliasResult(SUPAJNI.DDAPass_alias__SWIG_1(swigCPtr, this, SWIGTYPE_p_llvm__Value.getCPtr(V1), SWIGTYPE_p_llvm__Value.getCPtr(V2)), true);
  }

  public boolean runOnModule(SVFModule module) {
    return SUPAJNI.DDAPass_runOnModule__SWIG_0(swigCPtr, this, SVFModule.getCPtr(module), module);
  }

  public boolean runOnModule(SWIGTYPE_p_llvm__Module module) {
    return SUPAJNI.DDAPass_runOnModule__SWIG_1(swigCPtr, this, SWIGTYPE_p_llvm__Module.getCPtr(module));
  }

  public void selectClient(SVFModule module) {
    SUPAJNI.DDAPass_selectClient(swigCPtr, this, SVFModule.getCPtr(module), module);
  }

  public SWIGTYPE_p_llvm__StringRef getPassName() {
    return new SWIGTYPE_p_llvm__StringRef(SUPAJNI.DDAPass_getPassName(swigCPtr, this), true);
  }

  public void printQueryPTS() {
    SUPAJNI.DDAPass_printQueryPTS(swigCPtr, this);
  }

  public void runPointerAnalysis(SVFModule module, long kind) {
    SUPAJNI.DDAPass_runPointerAnalysis(swigCPtr, this, SVFModule.getCPtr(module), module, kind);
  }

  public void answerQueries(PointerAnalysis pta) {
    SUPAJNI.DDAPass_answerQueries(swigCPtr, this, PointerAnalysis.getCPtr(pta), pta);
  }

  public void initCxtInsensitiveEdges(PointerAnalysis pta, SWIGTYPE_p_SVFG svfg, SWIGTYPE_p_SCCDetectionT_SVFG_p_t svfgSCC, ConstSVFGEdgeSet insensitveEdges) {
    SUPAJNI.DDAPass_initCxtInsensitiveEdges(swigCPtr, this, PointerAnalysis.getCPtr(pta), pta, SWIGTYPE_p_SVFG.getCPtr(svfg), SWIGTYPE_p_SCCDetectionT_SVFG_p_t.getCPtr(svfgSCC), ConstSVFGEdgeSet.getCPtr(insensitveEdges), insensitveEdges);
  }

  public boolean edgeInSVFGSCC(SWIGTYPE_p_SCCDetectionT_SVFG_p_t svfgSCC, SWIGTYPE_p_SVFGEdge edge) {
    return SUPAJNI.DDAPass_edgeInSVFGSCC(swigCPtr, this, SWIGTYPE_p_SCCDetectionT_SVFG_p_t.getCPtr(svfgSCC), SWIGTYPE_p_SVFGEdge.getCPtr(edge));
  }

  public boolean edgeInCallGraphSCC(PointerAnalysis pta, SWIGTYPE_p_SVFGEdge edge) {
    return SUPAJNI.DDAPass_edgeInCallGraphSCC(swigCPtr, this, PointerAnalysis.getCPtr(pta), pta, SWIGTYPE_p_SVFGEdge.getCPtr(edge));
  }

  public void collectCxtInsenEdgeForRecur(PointerAnalysis pta, SWIGTYPE_p_SVFG svfg, ConstSVFGEdgeSet insensitveEdges) {
    SUPAJNI.DDAPass_collectCxtInsenEdgeForRecur(swigCPtr, this, PointerAnalysis.getCPtr(pta), pta, SWIGTYPE_p_SVFG.getCPtr(svfg), ConstSVFGEdgeSet.getCPtr(insensitveEdges), insensitveEdges);
  }

  public void collectCxtInsenEdgeForVFCycle(PointerAnalysis pta, SWIGTYPE_p_SVFG svfg, SWIGTYPE_p_SCCDetectionT_SVFG_p_t svfgSCC, ConstSVFGEdgeSet insensitveEdges) {
    SUPAJNI.DDAPass_collectCxtInsenEdgeForVFCycle(swigCPtr, this, PointerAnalysis.getCPtr(pta), pta, SWIGTYPE_p_SVFG.getCPtr(svfg), SWIGTYPE_p_SCCDetectionT_SVFG_p_t.getCPtr(svfgSCC), ConstSVFGEdgeSet.getCPtr(insensitveEdges), insensitveEdges);
  }

  public void set_pta(PointerAnalysis value) {
    SUPAJNI.DDAPass__pta_set(swigCPtr, this, PointerAnalysis.getCPtr(value), value);
  }

  public PointerAnalysis get_pta() {
    long cPtr = SUPAJNI.DDAPass__pta_get(swigCPtr, this);
    return (cPtr == 0) ? null : new PointerAnalysis(cPtr, false);
  }

  public void set_client(DDAClient value) {
    SUPAJNI.DDAPass__client_set(swigCPtr, this, DDAClient.getCPtr(value), value);
  }

  public DDAClient get_client() {
    long cPtr = SUPAJNI.DDAPass__client_get(swigCPtr, this);
    return (cPtr == 0) ? null : new DDAClient(cPtr, false);
  }

}
