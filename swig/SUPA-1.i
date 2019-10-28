%module SUPA
%{
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/FlowDDA.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAPass.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAVFSolver.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAStat.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/ContextDDA.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAClient.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/MemoryModel/PointerAnalysis.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/BasicTypes.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/PTAStat.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/SVFUtil.h"
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/SVFModule.h"
%}

%include <cpointer.i>
%include <carrays.i>
%include <stdint.i>
%include <stl.i>
%include <typemaps.i>
%include <std_set.i>
%include <std_map.i>
%include <std_string.i>
%include <std_pair.i>
%include <std_list.i>
%include <std_vector.i>
%include <std_deque.i>
#include <vector>
#include <stack>
#include <time.h>
#include <map>
#include <iostream>
#include <set>
#include <algorithm>
#include <deque>
#include <string>
#include <list>

%import <llvm/Pass.h>
#include <llvm/Bitcode/BitcodeReader.h>     /// for isBitcode
#include <llvm/Bitcode/BitcodeWriter.h>		// for WriteBitcodeToFile
#include <llvm/IR/CallSite.h>
#include <llvm/IR/InstIterator.h>	// for inst iteration
#include <llvm/IR/IRBuilder.h>		// for instrument svf.main
#include <llvm/IR/InstVisitor.h>	// for instruction visitor
#include <llvm/IR/GetElementPtrTypeIterator.h>	//for gep iterator
#include <llvm/IR/Instructions.h>
#include <llvm/IR/GlobalVariable.h>	// for GlobalVariable
#include <llvm/IR/DebugInfo.h>
#include <llvm/IRReader/IRReader.h>	// IR reader for bit file
#include <llvm/Support/GraphWriter.h>		// for graph write
#include <llvm/Support/SourceMgr.h> // for SMDiagnostic
#include <llvm/Support/raw_ostream.h>	// for output
#include <llvm/Transforms/Utils/Local.h>	// for FindDbgAddrUses
#include <llvm/Transforms/Utils/UnifyFunctionExitNodes.h>
#include <llvm/ADT/DenseSet.h>		// for dense map, set
#include <llvm/ADT/SmallVector.h>		// for small vector
#include <llvm/ADT/SparseBitVector.h>	// for points-to
#include <llvm/ADT/GraphTraits.h>		// for Graphtraits
#include <llvm/ADT/StringExtras.h>	// for utostr_32
#include <llvm/Analysis/ScalarEvolutionExpressions.h>
#include <llvm/Analysis/AliasAnalysis.h>
#include <llvm/Analysis/PostDominators.h>
#include <llvm/Analysis/DominanceFrontier.h>
#include <llvm/Analysis/CallGraph.h>	// call graph
#include <llvm/Analysis/ScalarEvolution.h>

#include "Util/CPPUtil.h"
#include "MSSA/SVFGBuilder.h"
#include "MemoryModel/PAGBuilder.h"
#include "WPA/CSC.h"
#include "MemoryModel/ConditionalPT.h"
#include "Util/GraphPrinter.h"
#include "Util/ThreadAPI.h"
#include "Util/ExtAPI.h"
#include "Util/DataFlowUtil.h"
#include "Util/PathCondAllocator.h"
#include "MemoryModel/PAG.h"
#include "Util/GEPTypeBridgeIterator.h"
#include "MemoryModel/PointsToDFDS.h"
#include "Util/DPItem.h"
#include "MSSA/SVFG.h"
#include "Util/PTACallGraph.h"
#include "WPA/Andersen.h"
#include "MemoryModel/PointsToDS.h"
#include "Util/Casting.h"
#include "Util/SCC.h"
//%pointer_class(PAG, PAGP)
///%pointer_class(PointerAnalysis, PointerAnalysisP)
//%pointer_class(Value, ValueP)
//%pointer_class(SVFGSCC, SVFGSCCP)
//%pointer_class(DDAClient, DDAClientP)
// %pointer_class(SVFGEdge, SVFGEdgeP)
// %pointer_class(SVFG, SVFGP)
// %pointer_class(PAGNode,PAGNodeP)
// %pointer_class(Function,FunctionP)
// %pointer_class(StoreSVFGNode,StoreSVFGNodeP)
// %pointer_class(MemObj,MemObjP)
// %pointer_class(AddrSVFGNode,AddrSVFGNodeP)
// %pointer_class(DDAStat,DDAStatP)
// %pointer_class(Instruction,InstructionP)
// %pointer_class(GepSVFGNode,GepSVFGNodeP)
// %pointer_class(AndersenWaveDiff,AndersenWaveDiffP)
// %pointer_class(PTACallGraph,PTACallGraphP)
// %pointer_class(DirectSVFGEdge,DirectSVFGEdgeP)
// %pointer_class(IndirectSVFGEdge,IndirectSVFGEdgeP)
// %pointer_class(CallGraphSCC,CallGraphSCCP)
// %pointer_class(LoadSVFGNode,LoadSVFGNodeP)
// %pointer_class(AllocaInst,AllocaInstP)
// %pointer_class(DiffPTDataTy,DiffPTDataTyP)
// %pointer_class(PTAStat,PTAStatP)
// %pointer_class(PTACallGraphNode,PTACallGraphNodeP)
// %pointer_class(ICFG,ICFGP)
// %pointer_class(PTDataTy,PTDataTyP)
// %pointer_class(CallInst,CallInstP)
// %pointer_class(IncDFPTDataTy,IncDFPTDataTyP)
// %pointer_class(CHGraph,CHGraphP)
// %pointer_class(TypeSystem,TypeSystemP)
// %pointer_class(LLVMModuleSet,LLVMModuleSetP)
// %pointer_class(LLVMContext,LLVMContextP)
// %pointer_class(Module,ModuleP)
// %pointer_class(GlobalAlias,GlobalAliasP)
// %pointer_class(GlobalVariable,GlobalVariableP)
// %pointer_class(u32_t,u32_tP)
// %pointer_class(ConstantExpr,ConstantExprP)
// %pointer_class(PointerType,PointerTypeP)
// %pointer_class(BasicBlock,BasicBlockP)
// %pointer_class(DominatorTree,DominatorTreeP)

namespace std {
%template(DPTItemSet) set<DPIm>;
%template(ConstSVFGEdgeSet) set<const SVFGEdge*>;
%template(DPImToCPtSetMap) map<DPIm, CPtSet>;
%template(DPMToCVarMap) map<DPIm,CVar>;
%template(DPMToDPMMap) map<DPIm,DPIm>;
%template(StoreToPMSetMap) map<const SVFGNode*, DPTItemSet>;
%template(CallSiteSet) set<CallSite>;
%template(FunctionSet) set<const Function*>;
%template(VTableSet) set<const GlobalValue*>;
%template(CallEdgeMap) map<CallSite, FunctionSet>;
%template(PtrToBVPtsMap) map<NodeID,PointsTo>;
%template(PtrCPtsMap) map<NodeID,CPtSet>;
%template(NodePair) pair<NodeID, NodeID>;
%template(NodeSet) set<NodeID>;
%template(NodeVector) vector<NodeID>;
%template(EdgeVector) vector<EdgeID>;
%template(NodeList) list<NodeID>;
%template(NodeDeque) deque<NodeID>;
%template(NUMStatMap) map<const char*,u32_t>;
%template(TIMEStatMap) map<const char*,double>;
%template(FunctionSetType) vector<Function*> ;
%template(GlobalSetType) vector<GlobalVariable*> ;
%template(AliasSetType) vector<GlobalAlias*> ;
%template(FunDeclToDefMapTy) map<const Function*, Function*>;
%template(FunDefToDeclsMapTy) map<const Function*, FunctionSetType> ;
%template(GlobalDefToRepMapTy) map<const GlobalVariable*, GlobalVariable*> ;
%template(StringVector) vector<std::string>;
%template(BasicBlockVector) vector<const BasicBlock*>;
%template(InstructionVector) vector<const Instruction*>;
}
typedef StmtDPItem<SVFGNode> LocDPItem;
typedef CxtStmtDPItem<SVFGNode> CxtLocDPItem;
typedef unsigned NodeID;
typedef unsigned EdgeID;
typedef unsigned SymID;
typedef unsigned CallSiteID;
typedef unsigned ThreadID;
typedef unsigned u32_t;
typedef unsigned long long u64_t;
typedef signed s32_t;
typedef signed long Size_t;
typedef llvm::SparseBitVector<> NodeBS;
typedef NodeBS PointsTo;
typedef PointsTo AliasSet;
typedef std::pair<NodeID, NodeID> NodePair;
typedef std::set<NodeID> NodeSet;
typedef llvm::DenseSet<NodePair,llvm::DenseMapInfo<std::pair<NodeID,NodeID> > > NodePairSet;
typedef llvm::DenseMap<NodePair,NodeID,llvm::DenseMapInfo<std::pair<NodeID,NodeID> > > NodePairMap;
typedef std::vector<NodeID> NodeVector;
typedef std::vector<EdgeID> EdgeVector;
typedef std::stack<NodeID> NodeStack;
typedef std::list<NodeID> NodeList;
typedef std::deque<NodeID> NodeDeque;
typedef llvm::SmallVector<u32_t,16> SmallVector16;
typedef llvm::SmallVector<u32_t,8> SmallVector8;
typedef NodeSet EdgeSet;
typedef SmallVector16 CallStrCxt;
typedef llvm::StringMap<u32_t> StringMap;
typedef llvm::SMDiagnostic SMDiagnostic;
typedef llvm::LLVMContext LLVMContext;
typedef llvm::Type Type;
typedef llvm::Function Function;
typedef llvm::BasicBlock BasicBlock;
typedef llvm::Value Value;
typedef llvm::Instruction Instruction;
typedef llvm::CallSite CallSite;
typedef llvm::GlobalValue GlobalValue;
typedef llvm::GlobalVariable GlobalVariable;
typedef llvm::Module Module;
typedef llvm::CallGraph LLVMCallGraph;
typedef llvm::User User;
typedef llvm::Use Use;
typedef llvm::Loop Loop;
typedef llvm::LoopInfo LoopInfo;
typedef llvm::UnifyFunctionExitNodes UnifyFunctionExitNodes;
typedef llvm::ModulePass ModulePass;
typedef llvm::AnalysisUsage AnalysisUsage;
typedef llvm::raw_ostream raw_ostream;
typedef llvm::raw_string_ostream raw_string_ostream;
typedef llvm::raw_fd_ostream raw_fd_ostream;
typedef llvm::StringRef StringRef;
typedef llvm::ToolOutputFile ToolOutputFile;
typedef llvm::StructType StructType;
typedef llvm::ArrayType ArrayType;
typedef llvm::PointerType PointerType;
typedef llvm::FunctionType FunctionType;
typedef llvm::VectorType VectorType;
typedef llvm::DataLayout DataLayout;
typedef llvm::StructLayout StructLayout;
typedef llvm::SmallVector<BasicBlock*, 8> SmallBBVector;
typedef llvm::ConstantStruct ConstantStruct;
typedef llvm::MemoryLocation MemoryLocation;
typedef llvm::Argument Argument;
typedef llvm::Constant Constant;
typedef llvm::ConstantData ConstantData;
typedef llvm::ConstantExpr ConstantExpr;
typedef llvm::ConstantAggregate ConstantAggregate;
typedef llvm::ConstantPointerNull ConstantPointerNull;
typedef llvm::ConstantArray ConstantArray;
typedef llvm::GlobalAlias GlobalAlias;
typedef llvm::AliasResult AliasResult;
typedef llvm::AnalysisID AnalysisID;
typedef llvm::ConstantDataArray ConstantDataArray;
typedef llvm::NamedMDNode NamedMDNode;
typedef llvm::MDString MDString;
typedef llvm::MDNode MDNode;
typedef llvm::AllocaInst AllocaInst;
typedef llvm::CallInst CallInst;
typedef llvm::InvokeInst InvokeInst;
typedef llvm::StoreInst StoreInst;
typedef llvm::LoadInst LoadInst;
typedef llvm::PHINode PHINode;
typedef llvm::GetElementPtrInst GetElementPtrInst;
typedef llvm::CastInst CastInst;
typedef llvm::BitCastInst BitCastInst;
typedef llvm::ReturnInst ReturnInst;
typedef llvm::ConstantInt ConstantInt;
typedef llvm::SelectInst SelectInst;
typedef llvm::IntToPtrInst IntToPtrInst;
typedef llvm::CmpInst CmpInst;
typedef llvm::BranchInst BranchInst;
typedef llvm::SwitchInst SwitchInst;
typedef llvm::ExtractValueInst  ExtractValueInst;
typedef llvm::InsertValueInst InsertValueInst;
typedef llvm::BinaryOperator BinaryOperator;
typedef llvm::PtrToIntInst PtrToIntInst;
typedef llvm::VAArgInst VAArgInst;
typedef llvm::ExtractElementInst ExtractElementInst;
typedef llvm::InsertElementInst InsertElementInst;
typedef llvm::ShuffleVectorInst ShuffleVectorInst;
typedef llvm::LandingPadInst LandingPadInst;
typedef llvm::ResumeInst ResumeInst;
typedef llvm::UnreachableInst UnreachableInst;
typedef llvm::FenceInst FenceInst;
typedef llvm::AtomicCmpXchgInst AtomicCmpXchgInst;
typedef llvm::AtomicRMWInst AtomicRMWInst;
typedef llvm::UndefValue UndefValue;
typedef llvm::FunctionCallee FunctionCallee;
typedef llvm::ScalarEvolutionWrapperPass ScalarEvolutionWrapperPass;
typedef llvm::ScalarEvolution ScalarEvolution;
typedef llvm::SCEVAddRecExpr SCEVAddRecExpr;
typedef llvm::SCEVConstant SCEVConstant;
typedef llvm::SCEV SCEV;
typedef llvm::DominanceFrontier DominanceFrontier;
typedef llvm::DominatorTree DominatorTree;
typedef llvm::PostDominatorTree PostDominatorTree;
typedef llvm::DomTreeNode DomTreeNode;
typedef llvm::DominanceFrontierBase<BasicBlock, false> DominanceFrontierBase;
typedef llvm::PostDominatorTreeWrapperPass PostDominatorTreeWrapperPass;
typedef llvm::inst_iterator inst_iterator;
typedef llvm::succ_const_iterator succ_const_iterator;
typedef llvm::const_inst_iterator const_inst_iterator;
typedef llvm::const_pred_iterator const_pred_iterator;
typedef llvm::gep_type_iterator gep_type_iterator;
typedef llvm::bridge_gep_iterator bridge_gep_iterator;
typedef llvm::GraphPrinter GraphPrinter;
typedef llvm::IRBuilder<> IRBuilder;
typedef llvm::IntegerType IntegerType;
class CHGraph;
class CHNode;

class TypeSystem;
class SVFModule;
class ICFG;
class PTAStat;
typedef DPItem DPIm;
/*
 * Pointer Analysis Base Class
 */
class PointerAnalysis {

public:
    /// Pointer analysis type list
    enum PTATY {
        // Whole program analysis
        Andersen_WPA,		///< Andersen PTA
        AndersenLCD_WPA,	///< Lazy cycle detection andersen-style WPA
        AndersenHCD_WPA,    ///< Hybird cycle detection andersen-style WPA
        AndersenHLCD_WPA,   ///< Hybird lazy cycle detection andersen-style WPA
        AndersenSCD_WPA,    ///< Selective cycle detection andersen-style WPA
        AndersenSFR_WPA,    ///< Stride-based field representation
        AndersenWaveDiff_WPA,	///< Diff wave propagation andersen-style WPA
        AndersenWaveDiffWithType_WPA,	///< Diff wave propagation with type info andersen-style WPA
        CSCallString_WPA,	///< Call string based context sensitive WPA
        CSSummary_WPA,		///< Summary based context sensitive WPA
        FSDATAFLOW_WPA,	///< Traditional Dataflow-based flow sensitive WPA
        FSSPARSE_WPA,		///< Sparse flow sensitive WPA
        FSCS_WPA,			///< Flow-, context- sensitive WPA
        FSCSPS_WPA,		///< Flow-, context-, path- sensitive WPA
        ADAPTFSCS_WPA,		///< Adaptive Flow-, context-, sensitive WPA
        ADAPTFSCSPS_WPA,	///< Adaptive Flow-, context-, path- sensitive WPA
        TypeCPP_WPA, ///<  Type-based analysis for C++

        // Demand driven analysis
        FieldS_DDA,		///< Field sensitive DDA
        FlowS_DDA,		///< Flow sensitive DDA
        PathS_DDA,		///< Guarded value-flow DDA
        Cxt_DDA,		///< context sensitive DDA


        Default_PTA		///< default pta without any analysis
    };

    /// Indirect call edges type, map a callsite to a set of callees
    //@{
    typedef llvm::AliasAnalysis AliasAnalysis;
    typedef std::set<CallSite> CallSiteSet;
    typedef PAG::CallSiteToFunPtrMap CallSiteToFunPtrMap;
    typedef	std::set<const Function*> FunctionSet;
    typedef std::map<CallSite, FunctionSet> CallEdgeMap;
    typedef SCCDetection<PTACallGraph*> CallGraphSCC;
    typedef std::set<const GlobalValue*> VTableSet;
    typedef std::set<const Function*> VFunSet;
    //@}

private:
    /// Release the memory
    void destroy();

protected:

    /// User input flags
    //@{
    /// Flag for printing the statistic results
    bool print_stat;
    /// Flag for iteration budget for on-the-fly statistics
    u32_t OnTheFlyIterBudgetForStat;
    //@}

    /// PAG
    static PAG* pag;
    /// Module
    SVFModule svfMod;
    /// Pointer analysis Type
    PTATY ptaTy;
    /// Statistics
    PTAStat* stat;
    /// Call graph used for pointer analysis
    PTACallGraph* ptaCallGraph;
    /// SCC for CallGraph
    CallGraphSCC* callGraphSCC;
    /// Interprocedural control-flow graph
    ICFG* icfg;
    /// CHGraph
    static CHGraph *chgraph;
    /// TypeSystem
    TypeSystem *typeSystem;

public:
    /// Return number of resolved indirect call edges
    inline Size_t getNumOfResolvedIndCallEdge() const {
        return getPTACallGraph()->getNumOfResolvedIndCallEdge();
    }
    /// Return call graph
    inline PTACallGraph* getPTACallGraph() const {
        return ptaCallGraph;
    }
    /// Return call graph SCC
    inline CallGraphSCC* getCallGraphSCC() const {
        return callGraphSCC;
    }

    /// Constructor
    PointerAnalysis(PTATY ty = Default_PTA);

    /// Type of pointer analysis
    inline PTATY getAnalysisTy() const {
        return ptaTy;
    }

    /// Get/set PAG
    ///@{
    inline PAG* getPAG() const {
        return pag;
    }
    static inline void setPAG(PAG* g) {
        pag = g;
    }
    //@}

    /// Get PTA stat
    inline PTAStat* getStat() const {
        return stat;
    }
    /// Module
    inline SVFModule getModule() const {
        return svfMod;
    }
    /// Get all Valid Pointers for resolution
    inline NodeSet& getAllValidPtrs() {
        return pag->getAllValidPtrs();
    }

    /// Destructor
    virtual ~PointerAnalysis();

    /// Initialization of a pointer analysis, including building symbol table and PAG etc.
    virtual void initialize(SVFModule svfModule);

    /// Finalization of a pointer analysis, including checking alias correctness
    virtual void finalize();

    /// Start Analysis here (main part of pointer analysis). It needs to be implemented in child class
    virtual void analyze(SVFModule svfModule) = 0;

    /// Compute points-to results on-demand, overridden by derived classes
    virtual void computeDDAPts(NodeID id) {}

    /// Interface exposed to users of our pointer analysis, given Location infos
    virtual AliasResult alias(const MemoryLocation &LocA,
                                    const MemoryLocation &LocB) = 0;

    /// Interface exposed to users of our pointer analysis, given Value infos
    virtual AliasResult alias(const Value* V1,
                                    const Value* V2) = 0;

    /// Interface exposed to users of our pointer analysis, given PAGNodeID
    virtual AliasResult alias(NodeID node1, NodeID node2) = 0;

protected:
    /// Return all indirect callsites
    inline const CallSiteToFunPtrMap& getIndirectCallsites() const {
        return pag->getIndirectCallsites();
    }
    /// Return function pointer PAGNode at a callsite cs
    inline NodeID getFunPtr(const CallSite& cs) const {
        return pag->getFunPtr(cs);
    }
    /// Alias check functions to verify correctness of pointer analysis
    //@{
    virtual void validateTests();
    virtual void validateSuccessTests(const char* fun);
    virtual void validateExpectedFailureTests(const char* fun);
    //@}

    /// Whether to dump the graph for debugging purpose
    bool dumpGraph();

    /// Reset all object node as field-sensitive.
    void resetObjFieldSensitive();

public:
    /// Dump the statistics
    void dumpStat();

    /// Determine whether a points-to contains a black hole or constant node
    //@{
    inline bool containBlackHoleNode(PointsTo& pts) {
        return pts.test(pag->getBlackHoleNode());
    }
    inline bool containConstantNode(PointsTo& pts) {
        return pts.test(pag->getConstantNode());
    }
    inline bool isBlkObjOrConstantObj(NodeID ptd) const {
        return pag->isBlkObjOrConstantObj(ptd);
    }
    inline bool isNonPointerObj(NodeID ptd) const {
        return pag->isNonPointerObj(ptd);
    }
    //@}

    /// Whether this object is heap or array
    //@{
    inline bool isHeapMemObj(NodeID id) const {
        const MemObj* mem = pag->getObject(id);
        assert(mem && "memory object is null??");
        return mem->isHeap();
    }

    inline bool isArrayMemObj(NodeID id) const {
        const MemObj* mem = pag->getObject(id);
        assert(mem && "memory object is null??");
        return mem->isArray();
    }
    //@}

    /// For field-sensitivity
    ///@{
    inline bool isFIObjNode(NodeID id) const {
        return (SVFUtil::isa<FIObjPN>(pag->getPAGNode(id)));
    }
    inline NodeID getBaseObjNode(NodeID id) {
        return pag->getBaseObjNode(id);
    }
    inline NodeID getFIObjNode(NodeID id) {
        return pag->getFIObjNode(id);
    }
    inline NodeID getGepObjNode(NodeID id, const LocationSet& ls) {
        return pag->getGepObjNode(id,ls);
    }
    inline const NodeBS& getAllFieldsObjNode(NodeID id) {
        return pag->getAllFieldsObjNode(id);
    }
    inline void setObjFieldInsensitive(NodeID id) {
        MemObj* mem =  const_cast<MemObj*>(pag->getBaseObj(id));
        mem->setFieldInsensitive();
    }
    inline bool isFieldInsensitive(NodeID id) const {
        const MemObj* mem =  pag->getBaseObj(id);
        return mem->isFieldInsensitive();
    }
    ///@}

    /// Whether print statistics
    inline bool printStat() {
        return print_stat;
    }

    /// Whether print statistics
    inline void disablePrintStat() {
        print_stat = false;
    }

    /// Get callees from an indirect callsite
    //@{
    inline CallEdgeMap& getIndCallMap() {
        return getPTACallGraph()->getIndCallMap();
    }
    inline bool hasIndCSCallees(CallSite cs) const {
        return getPTACallGraph()->hasIndCSCallees(cs);
    }
    inline const FunctionSet& getIndCSCallees(CallSite cs) const {
        return getPTACallGraph()->getIndCSCallees(cs);
    }
    inline const FunctionSet& getIndCSCallees(CallInst* csInst) const {
        CallSite cs = SVFUtil::getLLVMCallSite(csInst);
        return getIndCSCallees(cs);
    }
    //@}

    /// Resolve indirect call edges
    virtual void resolveIndCalls(CallSite cs, const PointsTo& target, CallEdgeMap& newEdges,LLVMCallGraph* callgraph = NULL);
    /// Match arguments for callsite at caller and callee
    inline bool matchArgs(CallSite cs, const Function* callee) {
        if(ThreadAPI::getThreadAPI()->isTDFork(cs))
            return true;
        else
            return cs.arg_size() == callee->arg_size();
    }

    /// CallGraph SCC related methods
    //@{
    /// CallGraph SCC detection
    inline void callGraphSCCDetection() {
        if(callGraphSCC==NULL)
            callGraphSCC = new CallGraphSCC(ptaCallGraph);

        callGraphSCC->find();
    }
    /// Get SCC rep node of a SVFG node.
    inline NodeID getCallGraphSCCRepNode(NodeID id) const {
        return callGraphSCC->repNode(id);
    }
    /// Return TRUE if this edge is inside a CallGraph SCC, i.e., src node and dst node are in the same SCC on the SVFG.
    inline bool inSameCallGraphSCC(const Function* fun1,const Function* fun2) {
        const PTACallGraphNode* src = ptaCallGraph->getCallGraphNode(fun1);
        const PTACallGraphNode* dst = ptaCallGraph->getCallGraphNode(fun2);
        return (getCallGraphSCCRepNode(src->getId()) == getCallGraphSCCRepNode(dst->getId()));
    }
    inline bool isInRecursion(const Function* fun) const {
        return callGraphSCC->isInCycle(ptaCallGraph->getCallGraphNode(fun)->getId());
    }
    /// Whether a local variable is in function recursions
    bool isLocalVarInRecursiveFun(NodeID id) const;
    //@}

    /// Return PTA name
    virtual const std::string PTAName() const {
        return "Pointer Analysis";
    }

    /// Get points-to targets of a pointer. It needs to be implemented in child class
    virtual PointsTo& getPts(NodeID ptr) = 0;
    
    /// Given an object, get all the nodes having whose pointsto contains the object. 
    /// Similar to getPts, this also needs to be implemented in child classes.
    virtual PointsTo& getRevPts(NodeID nodeId) = 0;

    /// Clear points-to data
    virtual void clearPts() {
    }

    /// Print targets of a function pointer
    void printIndCSTargets(const CallSite cs, const FunctionSet& targets);

    // Debug purpose
    //@{
    virtual void dumpTopLevelPtsTo() {}
    virtual void dumpAllPts() {}
    virtual void dumpCPts() {}
    virtual void dumpPts(NodeID ptr, const PointsTo& pts);
    void printIndCSTargets();
    void dumpAllTypes();
    //@}

    /// get CHGraph
    CHGraph *getCHGraph() const {
        return chgraph;
    }

    void getVFnsFromCHA(CallSite cs, std::set<const Function*> &vfns);
    void getVFnsFromPts(CallSite cs, const PointsTo &target, VFunSet &vfns);
    void connectVCallToVFns(CallSite cs, const VFunSet &vfns, CallEdgeMap& newEdges);
    virtual void resolveCPPIndCalls(CallSite cs,
                                    const PointsTo& target,
                                    CallEdgeMap& newEdges);

    /// get TypeSystem
    const TypeSystem *getTypeSystem() const {
        return typeSystem;
    }
};

/*!
 * Value-Flow Based Demand-Driven Points-to Analysis
 */
template<class CVar, class CPtSet, class DPIm>
class DDAVFSolver  {
    friend class DDAStat;
public:
    typedef SCCDetection<SVFG*> SVFGSCC;
    typedef SCCDetection<PTACallGraph*> CallGraphSCC;
    typedef PTACallGraphEdge::CallInstSet CallInstSet;
    typedef PAG::CallSiteSet CallSiteSet;
    typedef std::set<DPIm> DPTItemSet;
    typedef std::map<DPIm, CPtSet> DPImToCPtSetMap;
    typedef std::map<DPIm,CVar> DPMToCVarMap;
    typedef std::map<DPIm,DPIm> DPMToDPMMap;
    typedef llvm::DenseMap<NodeID, DPTItemSet> LocToDPMVecMap;
    typedef std::set<const SVFGEdge* > ConstSVFGEdgeSet;
    typedef SVFGEdge::SVFGEdgeSetTy SVFGEdgeSet;
    typedef std::map<const SVFGNode*, DPTItemSet> StoreToPMSetMap;

    ///Constructor
    DDAVFSolver(): outOfBudgetQuery(false),_pag(NULL),_svfg(NULL),_ander(NULL),_callGraph(NULL), _callGraphSCC(NULL), _svfgSCC(NULL), ddaStat(NULL) {
    }
    /// Destructor
    virtual ~DDAVFSolver() {
        if(_ander != NULL) {
            // AndersenWaveDiff::releaseAndersenWaveDiff();
            _ander = NULL;
        }

        if (_svfg != NULL) {
            // DDASVFGBuilder::releaseDDASVFG();
            _svfg = NULL;
        }

        if (_svfgSCC != NULL)
            delete _svfgSCC;
        _svfgSCC = NULL;

        _callGraph = NULL;
        _callGraphSCC = NULL;
    }
    /// Return candidate pointers for DDA
    inline NodeBS& getCandidateQueries() {
        return candidateQueries;
    }
    /// Given CVar and location (SVFGNode) return a new DPItem
    virtual inline DPIm getDPIm(const CVar& var, const SVFGNode* loc) const {
        DPIm dpm(var,loc);
        return dpm;
    }
    /// Union pts
    virtual bool unionDDAPts(CPtSet& pts, const CPtSet& targetPts) {
        return (pts |= targetPts);
    }
    /// Add pts
    virtual void addDDAPts(CPtSet& pts, const CVar& var) {
        pts.set(var);
    }
    /// Return SVFG
    inline SVFG* getSVFG() const {
        return _svfg;
    }
    /// Return SVFGSCC
    inline SVFGSCC* getSVFGSCC() const {
        return _svfgSCC;
    }
    // Dump cptsSet
    inline void dumpCPtSet(const CPtSet& cpts) const {
        SVFUtil::outs() << "{";
        for(typename CPtSet::iterator it = cpts.begin(), eit = cpts.end(); it!=eit; ++it) {
            SVFUtil::outs() << (*it) << " ";
        }
        SVFUtil::outs() << "}\n";
    }
    /// Compute points-to
    virtual const CPtSet& findPT(const DPIm& dpm) {

        if(isbkVisited(dpm)) {
            const CPtSet& cpts = getCachedPointsTo(dpm);
            DBOUT(DDDA, SVFUtil::outs() << "\t already backward visited dpm: ");
            DBOUT(DDDA, dpm.dump());
            DBOUT(DDDA, SVFUtil::outs() << "\t return points-to: ");
            DBOUT(DDDA, dumpCPtSet(cpts));
            return cpts;
        }

        DBOUT(DDDA, SVFUtil::outs() << "\t backward visit dpm: ");
        DBOUT(DDDA, dpm.dump());
        markbkVisited(dpm);
        addDpmToLoc(dpm);

        if(testOutOfBudget(dpm) == false) {

            CPtSet pts;
            handleSingleStatement(dpm, pts);

            /// Add successors of current stmt if its pts has been changed.
            updateCachedPointsTo(dpm, pts);
        }
        return getCachedPointsTo(dpm);
    }

protected:
    /// Handle single statement
    virtual void handleSingleStatement(const DPIm& dpm, CPtSet& pts) {
        /// resolve function pointer first at indirect callsite
        resolveFunPtr(dpm);

        const SVFGNode* node = dpm.getLoc();
        if(SVFUtil::isa<AddrSVFGNode>(node)) {
            handleAddr(pts,dpm,SVFUtil::cast<AddrSVFGNode>(node));
        }
        else if(SVFUtil::isa<CopySVFGNode>(node) || SVFUtil::isa<PHISVFGNode>(node)
                || SVFUtil::isa<ActualParmSVFGNode>(node) || SVFUtil::isa<FormalParmSVFGNode>(node)
                || SVFUtil::isa<ActualRetSVFGNode>(node) || SVFUtil::isa<FormalRetSVFGNode>(node)
                || SVFUtil::isa<NullPtrSVFGNode>(node)) {
            backtraceAlongDirectVF(pts,dpm);
        }
        else if(SVFUtil::isa<GepSVFGNode>(node)) {
            CPtSet gepPts;
            backtraceAlongDirectVF(gepPts,dpm);
            unionDDAPts(pts, processGepPts(SVFUtil::cast<GepSVFGNode>(node),gepPts));
        }
        else if(SVFUtil::isa<LoadSVFGNode>(node)) {
            CPtSet loadpts;
            startNewPTCompFromLoadSrc(loadpts,dpm);
            for(typename CPtSet::iterator it = loadpts.begin(), eit = loadpts.end(); it!=eit; ++it) {
                backtraceAlongIndirectVF(pts,getDPImWithOldCond(dpm,*it,node));
            }
        }
        else if(SVFUtil::isa<StoreSVFGNode>(node)) {
            if(isMustAlias(getLoadDpm(dpm),dpm)) {
                DBOUT(DDDA, SVFUtil::outs() << "+++must alias for load and store:");
                DBOUT(DDDA, getLoadDpm(dpm).dump());
                DBOUT(DDDA, dpm.dump());
                DBOUT(DDDA, SVFUtil::outs() << "+++\n");
                DOSTAT(ddaStat->_NumOfMustAliases++);
                backtraceToStoreSrc(pts,dpm);
            }
            else {
                CPtSet storepts;
                startNewPTCompFromStoreDst(storepts,dpm);
                for(typename CPtSet::iterator it = storepts.begin(), eit = storepts.end(); it!=eit; ++it) {
                    if(propagateViaObj(*it,getLoadCVar(dpm))) {
                        backtraceToStoreSrc(pts,getDPImWithOldCond(dpm,*it,node));

                        if(isStrongUpdate(storepts,SVFUtil::cast<StoreSVFGNode>(node))) {
                            DBOUT(DDDA, SVFUtil::outs() << "backward strong update for obj " << dpm.getCurNodeID() << "\n");
                            DOSTAT(addSUStat(dpm,node);)
                        }
                        else {
                            DOSTAT(rmSUStat(dpm,node);)
                            backtraceAlongIndirectVF(pts,getDPImWithOldCond(dpm,*it,node));
                        }
                    }
                    else {
                        backtraceAlongIndirectVF(pts,dpm);
                    }
                }
            }
        }
        else if(SVFUtil::isa<MRSVFGNode>(node)) {
            backtraceAlongIndirectVF(pts,dpm);
        }
        else
            assert(false && "unexpected kind of SVFG nodes");
    }

    /// recompute points-to for value-flow cycles and indirect calls
    void reCompute(const DPIm& dpm) {
        /// re-compute due to indirect calls
        SVFGEdgeSet newIndirectEdges;
        if(_pag->isFunPtr(dpm.getCurNodeID())) {
            const CallSiteSet& csSet = _pag->getIndCallSites(dpm.getCurNodeID());
            for(CallSiteSet::const_iterator it = csSet.begin(), eit = csSet.end(); it!=eit; ++it)
                updateCallGraphAndSVFG(dpm,*it,newIndirectEdges);
        }
        /// callgraph scc detection for local variable in recursion
        if(!newIndirectEdges.empty())
            _callGraphSCC->find();
        reComputeForEdges(dpm,newIndirectEdges,true);

        /// re-compute for transitive closures
        SVFGEdgeSet edgeSet(dpm.getLoc()->getOutEdges());
        reComputeForEdges(dpm,edgeSet,false);
    }

    /// Traverse along out edges to find all nodes which may be affected by locDPM.
    void reComputeForEdges(const DPIm& dpm, const SVFGEdgeSet& edgeSet, bool indirectCall = false) {
        for (SVFGNode::const_iterator it = edgeSet.begin(), eit = edgeSet.end(); it != eit; ++it) {
            const SVFGEdge* edge = *it;
            const SVFGNode* dst = edge->getDstNode();
            typename LocToDPMVecMap::const_iterator locIt = getLocToDPMVecMap().find(dst->getId());
            /// Only collect nodes we have traversed
            if (locIt == getLocToDPMVecMap().end())
                continue;
            DPTItemSet dpmSet(locIt->second.begin(), locIt->second.end());
            for(typename DPTItemSet::const_iterator it = dpmSet.begin(),eit = dpmSet.end(); it!=eit; ++it) {
                const DPIm& dstDpm = *it;
                if(!indirectCall && SVFUtil::isa<IndirectSVFGEdge>(edge) && !SVFUtil::isa<LoadSVFGNode>(edge->getDstNode())) {
                    if(dstDpm.getCurNodeID() == dpm.getCurNodeID()) {
                        DBOUT(DDDA,SVFUtil::outs() << "\t Recompute, forward from :");
                        DBOUT(DDDA, dpm.dump());
                        DOSTAT(ddaStat->_NumOfStepInCycle++);
                        clearbkVisited(dstDpm);
                        findPT(dstDpm);
                    }
                }
                else {
                    if(indirectCall)
                        DBOUT(DDDA,SVFUtil::outs() << "\t Recompute for indirect call from :");
                    else
                        DBOUT(DDDA,SVFUtil::outs() << "\t Recompute forward from :");
                    DBOUT(DDDA, dpm.dump());
                    DOSTAT(ddaStat->_NumOfStepInCycle++);
                    clearbkVisited(dstDpm);
                    findPT(dstDpm);
                }
            }
        }
    }

    /// Build SVFG
    virtual inline void buildSVFG(SVFModule module) {
        _ander = AndersenWaveDiff::createAndersenWaveDiff(module);
        _svfg = svfgBuilder.buildPTROnlySVFGWithoutOPT(_ander);
        _pag = _svfg->getPAG();
    }
    /// Reset visited map for next points-to query
    virtual inline void resetQuery() {
        if(outOfBudgetQuery)
            OOBResetVisited();

        locToDpmSetMap.clear();
        dpmToloadDpmMap.clear();
        loadToPTCVarMap.clear();
        outOfBudgetQuery = false;
        ddaStat->_NumOfStep = 0;
    }
    /// Reset visited map if the current query is out-of-budget
    inline void OOBResetVisited() {
        for(typename LocToDPMVecMap::const_iterator it = locToDpmSetMap.begin(),eit = locToDpmSetMap.end(); it!=eit; ++it) {
            DPTItemSet dpmSet(it->second.begin(), it->second.end());
            for(typename DPTItemSet::const_iterator dit = dpmSet.begin(),deit=dpmSet.end(); dit!=deit; ++dit)
                if(isOutOfBudgetDpm(*dit)==false)
                    clearbkVisited(*dit);
        }
    }
    /// GetDefinition SVFG
    inline const SVFGNode* getDefSVFGNode(const PAGNode* pagNode) const {
        return getSVFG()->getDefSVFGNode(pagNode);
    }
    /// Backward traverse along indirect value flows
    void backtraceAlongIndirectVF(CPtSet& pts, const DPIm& oldDpm) {
        const SVFGNode* node = oldDpm.getLoc();
        NodeID obj = oldDpm.getCurNodeID();
        if (_pag->isConstantObj(obj) || _pag->isNonPointerObj(obj))
            return;
        const SVFGEdgeSet edgeSet(node->getInEdges());
        for (SVFGNode::const_iterator it = edgeSet.begin(), eit = edgeSet.end(); it != eit; ++it) {
            if(const IndirectSVFGEdge* indirEdge = SVFUtil::dyn_cast<IndirectSVFGEdge>(*it)) {
                PointsTo& guard = const_cast<PointsTo&>(indirEdge->getPointsTo());
                if(guard.test(obj)) {
                    DBOUT(DDDA, SVFUtil::outs() << "\t\t==backtrace indirectVF svfgNode " <<
                          indirEdge->getDstID() << " --> " << indirEdge->getSrcID() << "\n");
                    backwardPropDpm(pts,oldDpm.getCurNodeID(),oldDpm,indirEdge);
                }
            }
        }
    }
    /// Backward traverse along direct value flows
    void backtraceAlongDirectVF(CPtSet& pts, const DPIm& oldDpm) {
        const SVFGNode* node = oldDpm.getLoc();
        const SVFGEdgeSet edgeSet(node->getInEdges());
        for (SVFGNode::const_iterator it = edgeSet.begin(), eit = edgeSet.end(); it != eit; ++it) {
            if(const DirectSVFGEdge* dirEdge = SVFUtil::dyn_cast<DirectSVFGEdge>(*it)) {
                DBOUT(DDDA, SVFUtil::outs() << "\t\t==backtrace directVF svfgNode " <<
                      dirEdge->getDstID() << " --> " << dirEdge->getSrcID() << "\n");
                const SVFGNode* srcNode = dirEdge->getSrcNode();
                backwardPropDpm(pts,getSVFG()->getLHSTopLevPtr(srcNode)->getId(),oldDpm,dirEdge);
            }
        }
    }

    /// Backward traverse for top-level pointers of load/store statements
    ///@{
    inline void startNewPTCompFromLoadSrc(CPtSet& pts, const DPIm& oldDpm) {
        const LoadSVFGNode* load = SVFUtil::cast<LoadSVFGNode>(oldDpm.getLoc());
        const SVFGNode* loadSrc = getDefSVFGNode(load->getPAGSrcNode());
        DBOUT(DDDA, SVFUtil::outs() << "!##start new computation from loadSrc svfgNode " <<
              load->getId() << " --> " << loadSrc->getId() << "\n");
        const SVFGEdge* edge = getSVFG()->getSVFGEdge(loadSrc,load,SVFGEdge::IntraDirectVF);
        assert(edge && "Edge not found!!");
        backwardPropDpm(pts,load->getPAGSrcNodeID(),oldDpm,edge);

    }
    inline void startNewPTCompFromStoreDst(CPtSet& pts, const DPIm& oldDpm) {
        const StoreSVFGNode* store = SVFUtil::cast<StoreSVFGNode>(oldDpm.getLoc());
        const SVFGNode* storeDst = getDefSVFGNode(store->getPAGDstNode());
        DBOUT(DDDA, SVFUtil::outs() << "!##start new computation from storeDst svfgNode " <<
              store->getId() << " --> " << storeDst->getId() << "\n");
        const SVFGEdge* edge = getSVFG()->getSVFGEdge(storeDst,store,SVFGEdge::IntraDirectVF);
        assert(edge && "Edge not found!!");
        backwardPropDpm(pts,store->getPAGDstNodeID(),oldDpm,edge);
    }
    inline void backtraceToStoreSrc(CPtSet& pts, const DPIm& oldDpm) {
        const StoreSVFGNode* store = SVFUtil::cast<StoreSVFGNode>(oldDpm.getLoc());
        const SVFGNode* storeSrc = getDefSVFGNode(store->getPAGSrcNode());
        DBOUT(DDDA, SVFUtil::outs() << "++backtrace to storeSrc from svfgNode " << getLoadDpm(oldDpm).getLoc()->getId() << " to "<<
              store->getId() << " to " << storeSrc->getId() <<"\n");
        const SVFGEdge* edge = getSVFG()->getSVFGEdge(storeSrc,store,SVFGEdge::IntraDirectVF);
        assert(edge && "Edge not found!!");
        backwardPropDpm(pts,store->getPAGSrcNodeID(),oldDpm,edge);
    }
    //@}

    /// dpm transit during backward tracing
    virtual void backwardPropDpm(CPtSet& pts, NodeID ptr,const DPIm& oldDpm,const SVFGEdge* edge) {
        DPIm dpm(oldDpm);
        dpm.setLocVar(edge->getSrcNode(),ptr);
        DOTIMESTAT(double start = DDAStat::getClk());
        /// handle context-/path- sensitivity
        if(handleBKCondition(dpm,edge)==false) {
            DOTIMESTAT(ddaStat->_TotalTimeOfBKCondition += DDAStat::getClk() - start);
            DBOUT(DDDA, SVFUtil::outs() << "\t!!! infeasible path svfgNode: " << edge->getDstID() << " --| " << edge->getSrcID() << "\n");
            DOSTAT(ddaStat->_NumOfInfeasiblePath++);
            return;
        }

        /// record the source of load dpm
        if(SVFUtil::isa<IndirectSVFGEdge>(edge))
            addLoadDpmAndCVar(dpm,getLoadDpm(oldDpm),getLoadCVar(oldDpm));

        DOSTAT(ddaStat->_NumOfDPM++);
        /// handle out of budget case
        unionDDAPts(pts,findPT(dpm));
    }
    /// whether load and store are aliased
    virtual bool isMustAlias(const DPIm& loadDpm, const DPIm& storeDPm) {
        return false;
    }
    /// Return TRUE if this is a strong update STORE statement.
    virtual bool isStrongUpdate(const CPtSet& dstCPSet, const StoreSVFGNode* store) {
        if (dstCPSet.count() == 1) {
            /// Find the unique element in cpts
            typename CPtSet::iterator it = dstCPSet.begin();
            const CVar& var = *it;
            // Strong update can be made if this points-to target is not heap, array or field-insensitive.
            if (!isHeapCondMemObj(var,store) && !isArrayCondMemObj(var)
                    && !isFieldInsenCondMemObj(var) && !isLocalCVarInRecursion(var)) {
                return true;
            }
        }
        return false;
    }
    /// Whether a local variable is in function recursions
    virtual inline bool isLocalCVarInRecursion(const CVar& var) const {
        NodeID id = getPtrNodeID(var);
        const MemObj* obj = _pag->getObject(id);
        assert(obj && "object not found!!");
        if(obj->isStack()) {
            if(const AllocaInst* local = SVFUtil::dyn_cast<AllocaInst>(obj->getRefVal())) {
                const Function* fun = local->getParent()->getParent();
                return _callGraphSCC->isInCycle(_callGraph->getCallGraphNode(fun)->getId());
            }
        }
        return false;
    }

    /// If the points-to contain the object obj, we could move forward along indirect value-flow edge
    virtual inline bool propagateViaObj(const CVar& storeObj, const CVar& loadObj) {
        if(getPtrNodeID(storeObj) == getPtrNodeID(loadObj))
            return true;
        return false;
    }
    /// resolve function pointer
    void resolveFunPtr(const DPIm& dpm) {
        if(Instruction* callInst= getSVFG()->isCallSiteRetSVFGNode(dpm.getLoc())) {
            CallSite cs = SVFUtil::getLLVMCallSite(callInst);
            if(_pag->isIndirectCallSites(cs)) {
                NodeID funPtr = _pag->getFunPtr(cs);
                DPIm funPtrDpm(dpm);
                funPtrDpm.setLocVar(getSVFG()->getDefSVFGNode(_pag->getPAGNode(funPtr)),funPtr);
                findPT(funPtrDpm);
            }
        }
        else if(const Function* fun = getSVFG()->isFunEntrySVFGNode(dpm.getLoc())) {
            CallInstSet csSet;
            /// use pre-analysis call graph to approximate all potential callsites
            _ander->getPTACallGraph()->getIndCallSitesInvokingCallee(fun,csSet);
            for(CallInstSet::const_iterator it = csSet.begin(), eit = csSet.end(); it!=eit; ++it) {
                CallSite cs = SVFUtil::getLLVMCallSite(*it);
                NodeID funPtr = _pag->getFunPtr(cs);
                DPIm funPtrDpm(dpm);
                funPtrDpm.setLocVar(getSVFG()->getDefSVFGNode(_pag->getPAGNode(funPtr)),funPtr);
                findPT(funPtrDpm);
            }
        }
    }
    /// Methods to be implemented in child class
    //@{
    /// Get variable ID (PAGNodeID) according to CVar
    virtual NodeID getPtrNodeID(const CVar& var) const = 0;
    /// ProcessGep node to generate field object nodes of a struct
    virtual CPtSet processGepPts(const GepSVFGNode* gep, const CPtSet& srcPts) = 0;
    /// Handle AddrSVFGNode to add proper points-to
    virtual void handleAddr(CPtSet& pts,const DPIm& dpm,const AddrSVFGNode* addr) = 0;
    /// Get conservative points-to results when the query is out of budget
    virtual CPtSet getConservativeCPts(const DPIm& dpm) = 0;
    /// Handle condition for context or path analysis (backward analysis)
    virtual inline bool handleBKCondition(DPIm& oldDpm, const SVFGEdge* edge) {
        return true;
    }
    /// Update call graph
    virtual inline void updateCallGraphAndSVFG(const DPIm& dpm,CallSite cs,SVFGEdgeSet& svfgEdges) {}
    //@}

    ///Visited flags to avoid cycles
    //@{
    inline void markbkVisited(const DPIm& dpm) {
        backwardVisited.insert(dpm);
    }
    inline bool isbkVisited(const DPIm& dpm) {
        return backwardVisited.find(dpm)!=backwardVisited.end();
    }
    inline void clearbkVisited(const DPIm& dpm) {
        assert(backwardVisited.find(dpm)!=backwardVisited.end() && "dpm not found!");
        backwardVisited.erase(dpm);
    }
    //@}

    /// Points-to Caching for top-level pointers and address-taken objects
    //@{
    virtual inline CPtSet& getCachedPointsTo(const DPIm& dpm) {
        if (isTopLevelPtrStmt(dpm.getLoc()))
            return getCachedTLPointsTo(dpm);
        else
            return getCachedADPointsTo(dpm);
    }
    virtual inline void updateCachedPointsTo(const DPIm& dpm, CPtSet& pts) {
        CPtSet& dpmPts = getCachedPointsTo(dpm);
        if (unionDDAPts(dpmPts, pts)) {
            DOSTAT(double start = DDAStat::getClk());
            reCompute(dpm);
            DOSTAT(ddaStat->_AnaTimeCyclePerQuery += DDAStat::getClk() - start);
        }
    }
    virtual inline CPtSet& getCachedTLPointsTo(const DPIm& dpm) {
        return dpmToTLCPtSetMap[dpm];
    }
    virtual inline CPtSet& getCachedADPointsTo(const DPIm& dpm) {
        return dpmToADCPtSetMap[dpm];
    }
    //@}

    /// Whether this is a top-level pointer statement
    inline bool isTopLevelPtrStmt(const SVFGNode* stmt) {
        if (SVFUtil::isa<StoreSVFGNode>(stmt) || SVFUtil::isa<MRSVFGNode>(stmt))
            return false;
        else
            return true;
    }
    /// Return dpm with old context and path conditions
    virtual inline DPIm getDPImWithOldCond(const DPIm& oldDpm,const CVar& var, const SVFGNode* loc) {
        DPIm dpm(oldDpm);
        dpm.setLocVar(loc,getPtrNodeID(var));

        if(SVFUtil::isa<StoreSVFGNode>(loc))
            addLoadDpmAndCVar(dpm,getLoadDpm(oldDpm),var);

        if(SVFUtil::isa<LoadSVFGNode>(loc))
            addLoadDpmAndCVar(dpm,oldDpm,var);

        DOSTAT(ddaStat->_NumOfDPM++);
        return dpm;
    }
    /// SVFG SCC detection
    inline void SVFGSCCDetection() {
        if(_svfgSCC==NULL) {
            _svfgSCC = new SVFGSCC(getSVFG());
        }
        _svfgSCC->find();
    }
    /// Get SCC rep node of a SVFG node.
    inline NodeID getSVFGSCCRepNode(NodeID id) {
        return _svfgSCC->repNode(id);
    }
    /// Return whether this SVFGNode is in cycle
    inline bool isSVFGNodeInCycle(const SVFGNode* node) {
        return _svfgSCC->isInCycle(node->getId());
    }
    /// Return TRUE if this edge is inside a SVFG SCC, i.e., src node and dst node are in the same SCC on the SVFG.
    inline bool edgeInSVFGSCC(const SVFGEdge* edge) {
        return (getSVFGSCCRepNode(edge->getSrcID()) == getSVFGSCCRepNode(edge->getDstID()));
    }
    /// Set callgraph
    inline void setCallGraph (PTACallGraph* cg) {
        _callGraph = cg;
    }
    /// Set callgraphSCC
    inline void setCallGraphSCC (CallGraphSCC* scc) {
        _callGraphSCC = scc;
    }
    /// Check heap and array object
    //@{
    virtual inline bool isHeapCondMemObj(const CVar& var, const StoreSVFGNode* store) {
        const MemObj* mem = _pag->getObject(getPtrNodeID(var));
        assert(mem && "memory object is null??");
        return mem->isHeap();
    }

    inline bool isArrayCondMemObj(const CVar& var) const {
        const MemObj* mem = _pag->getObject(getPtrNodeID(var));
        assert(mem && "memory object is null??");
        return mem->isArray();
    }
    inline bool isFieldInsenCondMemObj(const CVar& var) const {
        const MemObj* mem =  _pag->getBaseObj(getPtrNodeID(var));
        return mem->isFieldInsensitive();
    }
    //@}
private:
    /// Map a SVFGNode to its dpms for handling value-flow cycles
    //@{
    inline const LocToDPMVecMap& getLocToDPMVecMap() const {
        return locToDpmSetMap;
    }
    inline const DPTItemSet& getDpmSetAtLoc(const SVFGNode* loc) {
        return locToDpmSetMap[loc->getId()];
    }
    inline void addDpmToLoc(const DPIm& dpm) {
        locToDpmSetMap[dpm.getLoc()->getId()].insert(dpm);
    }
    inline void removeDpmFromLoc(const DPIm& dpm) {
        assert(dpm == locToDpmSetMap[dpm.getLoc()].back() && "dpm not match with the end of vector");
        locToDpmSetMap[dpm.getLoc()->getId()].erase(dpm);
    }
    //@}
protected:
    /// LoadDpm for must-alias analysis
    //@{
    inline void addLoadDpmAndCVar(const DPIm& dpm,const DPIm& loadDpm,const CVar& loadVar) {
        addLoadCVar(dpm,loadVar);
        addLoadDpm(dpm,loadDpm);
    }
    /// Note that simply use "dpmToloadDpmMap[dpm]=loadDpm", requires DPIm have a default constructor
    inline void addLoadDpm(const DPIm& dpm,const DPIm& loadDpm) {
        typename DPMToDPMMap::iterator it = dpmToloadDpmMap.find(dpm);
        if(it!=dpmToloadDpmMap.end())
            it->second = loadDpm;
        else
            dpmToloadDpmMap.insert(std::make_pair(dpm,loadDpm));
    }
    inline const DPIm& getLoadDpm(const DPIm& dpm) const {
        typename DPMToDPMMap::const_iterator it = dpmToloadDpmMap.find(dpm);
        assert(it!=dpmToloadDpmMap.end() && "not found??");
        return it->second;
    }
    inline void addLoadCVar(const DPIm& dpm, const CVar& loadVar) {
        typename DPMToCVarMap::iterator it = loadToPTCVarMap.find(dpm);
        if(it!=loadToPTCVarMap.end())
            it->second = loadVar;
        else
            loadToPTCVarMap.insert(std::make_pair(dpm,loadVar));
    }
    inline const CVar& getLoadCVar(const DPIm& dpm) const {
        typename DPMToCVarMap::const_iterator it = loadToPTCVarMap.find(dpm);
        assert(it!=loadToPTCVarMap.end() && "not found??");
        return it->second;
    }
    //@}
    /// Return Andersen's analysis
    inline AndersenWaveDiff* getAndersenAnalysis() const {
        return _ander;
    }
    /// handle out-of-budget queries
    //@{
    /// Handle out-of-budget dpm
    inline void handleOutOfBudgetDpm(const DPIm& dpm) {}
    inline bool testOutOfBudget(const DPIm& dpm) {
        if(outOfBudgetQuery) return true;
        if(++ddaStat->_NumOfStep > DPIm::getMaxBudget())
            outOfBudgetQuery = true;
        return isOutOfBudgetDpm(dpm) || outOfBudgetQuery;
    }
    inline bool isOutOfBudgetQuery() const {
        return outOfBudgetQuery;
    }
    inline void addOutOfBudgetDpm(const DPIm& dpm) {
        outOfBudgetDpms.insert(dpm);
    }
    inline bool isOutOfBudgetDpm(const DPIm& dpm) const {
        return outOfBudgetDpms.find(dpm) != outOfBudgetDpms.end();
    }
    //@}

    /// Set DDAStat
    inline DDAStat* setDDAStat(DDAStat* s) {
        ddaStat = s;
        return ddaStat;
    }
    /// stat strong updates num
    inline void addSUStat(const DPIm& dpm, const SVFGNode* node) {
        if (storeToDPMs[node].insert(dpm).second) {
            ddaStat->_NumOfStrongUpdates++;
            ddaStat->_StrongUpdateStores.set(node->getId());
        }
    }
    /// remove strong updates num if the dpm goes to weak updates branch
    inline void rmSUStat(const DPIm& dpm, const SVFGNode* node) {
        DPTItemSet& dpmSet = storeToDPMs[node];
        if (dpmSet.erase(dpm)) {
            ddaStat->_NumOfStrongUpdates--;
            if(dpmSet.empty())
                ddaStat->_StrongUpdateStores.reset(node->getId());
        }
    }

    bool outOfBudgetQuery;			///< Whether the current query is out of step limits
    PAG* _pag;						///< PAG
    SVFG* _svfg;					///< SVFG
    AndersenWaveDiff* _ander;		///< Andersen's analysis
    NodeBS candidateQueries;		///< candidate pointers;
    PTACallGraph* _callGraph;		///< CallGraph
    CallGraphSCC* _callGraphSCC;	///< SCC for CallGraph
    SVFGSCC* _svfgSCC;				///< SCC for SVFG
    DPTItemSet backwardVisited;		///< visited map during backward traversing
    DPImToCPtSetMap dpmToTLCPtSetMap;	///< points-to caching map for top-level vars
    DPImToCPtSetMap dpmToADCPtSetMap;	///< points-to caching map for address-taken vars
    LocToDPMVecMap locToDpmSetMap;	///< map location to its dpms
    DPMToDPMMap dpmToloadDpmMap;		///< dpms at loads for may/must-alias analysis with stores
    DPMToCVarMap loadToPTCVarMap;	///< map a load dpm to its cvar pointed by its pointer operand
    DPTItemSet outOfBudgetDpms;		///< out of budget dpm set
    StoreToPMSetMap storeToDPMs;	///< map store to set of DPM which have been stong updated there
    DDAStat* ddaStat;				///< DDA stat
    SVFGBuilder svfgBuilder;			///< SVFG Builder
};

/*!
 * Pointer analysis implementation which uses bit vector based points-to data structure
 */
class BVDataPTAImpl : public PointerAnalysis {

public:
    typedef PTData<NodeID,PointsTo> PTDataTy;	/// Points-to data structure type
    typedef DiffPTData<NodeID,PointsTo,EdgeID> DiffPTDataTy;	/// Points-to data structure type
    typedef DFPTData<NodeID,PointsTo> DFPTDataTy;	/// Points-to data structure type
    typedef IncDFPTData<NodeID,PointsTo> IncDFPTDataTy;	/// Points-to data structure type

    /// Constructor
    BVDataPTAImpl(PointerAnalysis::PTATY type);

    /// Destructor
    virtual ~BVDataPTAImpl() {
        destroy();
    }

    /// Release memory
    inline void destroy() {
        delete ptD;
        ptD = NULL;
    }

    /// Get points-to and reverse points-to
    ///@{
    virtual inline PointsTo& getPts(NodeID id) {
        return ptD->getPts(id);
    }
    virtual inline PointsTo& getRevPts(NodeID nodeId) {
        return ptD->getRevPts(nodeId);
    }
    //@}

    /// Expand FI objects
    void expandFIObjs(const PointsTo& pts, PointsTo& expandedPts);

    /// Interface for analysis result storage on filesystem.
    //@{
    virtual void writeToFile(const std::string& filename);
    virtual bool readFromFile(const std::string& filename);
    //@}

protected:

    /// Update callgraph. This should be implemented by its subclass.
    virtual inline bool updateCallGraph(const CallSiteToFunPtrMap& callsites) {
        assert(false && "Virtual function not implemented!");
        return false;
    }

    /// Get points-to data structure
    inline PTDataTy* getPTDataTy() const {
        return ptD;
    }
    inline DiffPTDataTy* getDiffPTDataTy() const {
        return SVFUtil::cast<DiffPTDataTy>(ptD);
    }
    inline IncDFPTDataTy* getDFPTDataTy() const {
        return SVFUtil::cast<IncDFPTDataTy>(ptD);
    }

    /// Union/add points-to. Add the reverse points-to for node collapse purpose
    /// To be noted that adding reverse pts might incur 10% total overhead during solving
    //@{
    virtual inline bool unionPts(NodeID id, const PointsTo& target) {
        return ptD->unionPts(id, target);
    }
    virtual inline bool unionPts(NodeID id, NodeID ptd) {
        return ptD->unionPts(id,ptd);
    }
    virtual inline bool addPts(NodeID id, NodeID ptd) {
        return ptD->addPts(id,ptd);
    }
    //@}

    /// Clear all data
    virtual inline void clearPts() {
        ptD->clear();
    }

    /// On the fly call graph construction
    virtual void onTheFlyCallGraphSolve(const CallSiteToFunPtrMap& callsites, CallEdgeMap& newEdges);

private:
    /// Points-to data
    PTDataTy* ptD;

public:
    /// Interface expose to users of our pointer analysis, given Location infos
    virtual AliasResult alias(const MemoryLocation  &LocA,
                                    const MemoryLocation  &LocB);

    /// Interface expose to users of our pointer analysis, given Value infos
    virtual AliasResult alias(const Value* V1,
                                    const Value* V2);

    /// Interface expose to users of our pointer analysis, given PAGNodeID
    virtual AliasResult alias(NodeID node1, NodeID node2);

    /// Interface expose to users of our pointer analysis, given two pts
    virtual AliasResult alias(const PointsTo& pts1, const PointsTo& pts2);

    /// dump and debug, print out conditional pts
    //@{
    virtual void dumpCPts() {
        ptD->dumpPTData();
    }

    virtual void dumpTopLevelPtsTo();

    virtual void dumpAllPts();
    //@}
};
class DDAClient;

/*!
 * Flow sensitive demand-driven analysis on value-flow graph
 */
%template(DDAVFSolverNPL) DDAVFSolver<NodeID,PointsTo,LocDPItem>;
class FlowDDA : public BVDataPTAImpl, public DDAVFSolver<NodeID,PointsTo,LocDPItem> {

public:
    typedef BVDataPTAImpl::CallSiteSet CallSiteSet;
    typedef BVDataPTAImpl::CallEdgeMap	CallEdgeMap;
    typedef BVDataPTAImpl::FunctionSet	FunctionSet;
    /// Constructor
    FlowDDA(SVFModule m, DDAClient* client): BVDataPTAImpl(PointerAnalysis::FlowS_DDA),
        DDAVFSolver<NodeID,PointsTo,LocDPItem>(),
        _client(client) {
    }
    /// Destructor
    inline virtual ~FlowDDA() {
    }
    /// dummy analyze method
    virtual void analyze(SVFModule mod) {}

    /// Compute points-to set for all top variable
    void computeDDAPts(NodeID id);

    /// Handle out-of-budget dpm
    void handleOutOfBudgetDpm(const LocDPItem& dpm);

    /// Handle condition for flow analysis (backward analysis)
    virtual bool handleBKCondition(LocDPItem& dpm, const SVFGEdge* edge);

    /// refine indirect call edge
    bool testIndCallReachability(LocDPItem& dpm, const Function* callee, CallSiteID csId);

    /// Initialization of the analysis
    inline virtual void initialize(SVFModule module) {
        BVDataPTAImpl::initialize(module);
        buildSVFG(module);
        setCallGraph(getPTACallGraph());
        setCallGraphSCC(getCallGraphSCC());
        stat = setDDAStat(new DDAStat(this));
    }

    /// Finalize analysis
    inline virtual void finalize() {
        BVDataPTAImpl::finalize();
    }

    /// we exclude concrete heap here following the conditions:
    /// (1) local allocated heap and
    /// (2) not escaped to the scope outside the current function
    /// (3) not inside loop
    /// (4) not involved in recursion
    bool isHeapCondMemObj(const NodeID& var, const StoreSVFGNode* store);

    /// Override parent method
    inline PointsTo getConservativeCPts(const LocDPItem& dpm) {
        return getAndersenAnalysis()->getPts(dpm.getCurNodeID());
    }
    /// Override parent method
    virtual inline NodeID getPtrNodeID(const NodeID& var) const {
        return var;
    }
    /// Handle Address SVFGNode to add proper points-to
    inline void handleAddr(PointsTo& pts,const LocDPItem& dpm,const AddrSVFGNode* addr) {
        NodeID srcID = addr->getPAGSrcNodeID();
        /// whether this object is set field-insensitive during pre-analysis
        if (isFieldInsensitive(srcID))
            srcID = getFIObjNode(srcID);

        addDDAPts(pts,srcID);
        DBOUT(DDDA, SVFUtil::outs() << "\t add points-to target " << srcID << " to dpm ");
        DBOUT(DDDA, dpm.dump());
    }
    /// processGep node
    PointsTo processGepPts(const GepSVFGNode* gep, const PointsTo& srcPts);

    /// Update call graph.
    //@{
    void updateCallGraphAndSVFG(const LocDPItem& dpm,CallSite cs,SVFGEdgeSet& svfgEdges)
    {
        CallEdgeMap newEdges;
        resolveIndCalls(cs, getCachedPointsTo(dpm), newEdges);
        for (CallEdgeMap::const_iterator iter = newEdges.begin(),eiter = newEdges.end(); iter != eiter; iter++) {
            CallSite newcs = iter->first;
            const FunctionSet & functions = iter->second;
            for (FunctionSet::const_iterator func_iter = functions.begin(); func_iter != functions.end(); func_iter++) {
                const Function * func = *func_iter;
                getSVFG()->connectCallerAndCallee(newcs, func, svfgEdges);
            }
        }
    }
    //@}

    virtual const std::string PTAName() const {
        return "FlowSensitive DDA";
    }

private:

    /// Override parent class functions to get/add cached points-to directly via PAGNode ID
    //@{
    inline PointsTo& getCachedTLPointsTo(const LocDPItem& dpm) {
        return getPts(dpm.getCurNodeID());
    }
    //@}

    DDAClient* _client;				///< DDA client
    PTACFInfoBuilder loopInfoBuilder; ///< LoopInfo
};

/*
 * @file: DDAPass.h
 * @author: Yulei Sui
 * @date: 01/07/2014
 * @version: 1.0
 *
 */

/*!
 * Demand-Driven Pointer Analysis.
 * This class performs various pointer analysis on the given module.
 */
class DDAPass: public ModulePass {

public:
    /// Pass ID
    static char ID;
    typedef SCCDetection<SVFG*> SVFGSCC;
    typedef std::set<const SVFGEdge*> SVFGEdgeSet;
    typedef std::vector<PointerAnalysis*> PTAVector;

    DDAPass() : ModulePass(ID), _pta(NULL), _client(NULL) {}
    ~DDAPass();

    virtual inline void getAnalysisUsage(AnalysisUsage &au) const {
        // declare your dependencies here.
        /// do not intend to change the IR in this pass,
        au.setPreservesAll();
    }

    virtual inline void* getAdjustedAnalysisPointer(AnalysisID id) {
        return this;
    }

    /// Interface expose to users of our pointer analysis, given Location infos
    virtual inline AliasResult alias(const MemoryLocation &LocA, const MemoryLocation &LocB) {
        return alias(LocA.Ptr, LocB.Ptr);
    }

    /// Interface expose to users of our pointer analysis, given Value infos
    virtual AliasResult alias(const Value* V1,	const Value* V2);

    /// We start from here
    virtual bool runOnModule(SVFModule module);

    /// We start from here
    virtual bool runOnModule(Module& module) {
        return runOnModule(module);
    }

    /// Select a client
    virtual void selectClient(SVFModule module);

    /// Pass name
    virtual inline StringRef getPassName() const {
        return "DDAPass";
    }

    /// Print queries' pts
    void printQueryPTS();
    /// Create pointer analysis according to specified kind and analyze the module.
    void runPointerAnalysis(SVFModule module, u32_t kind);
    /// Initialize queries for DDA
    void answerQueries(PointerAnalysis* pta);
    /// Context insensitive Edge for DDA
    void initCxtInsensitiveEdges(PointerAnalysis* pta, const SVFG* svfg,const SVFGSCC* svfgSCC, SVFGEdgeSet& insensitveEdges);
    /// Return TRUE if this edge is inside a SVFG SCC, i.e., src node and dst node are in the same SCC on the SVFG.
    bool edgeInSVFGSCC(const SVFGSCC* svfgSCC,const SVFGEdge* edge);
    /// Return TRUE if this edge is inside a SVFG SCC, i.e., src node and dst node are in the same SCC on the SVFG.
    bool edgeInCallGraphSCC(PointerAnalysis* pta,const SVFGEdge* edge);

    void collectCxtInsenEdgeForRecur(PointerAnalysis* pta, const SVFG* svfg,SVFGEdgeSet& insensitveEdges);
    void collectCxtInsenEdgeForVFCycle(PointerAnalysis* pta, const SVFG* svfg,const SVFGSCC* svfgSCC, SVFGEdgeSet& insensitveEdges);

    PointerAnalysis* _pta;	///<  pointer analysis to be executed.
    DDAClient* _client;		///<  DDA client used

};

class PointerAnalysis;

/*!
 * Pointer Analysis Statistics
 */
class PTAStat {
public:
    static const char* TotalAnalysisTime; ///< Total analysis time
    static const char* SCCDetectionTime; ///< Total SCC detection time
    static const char* SCCMergeTime; ///< Total SCC merge time

    static const char* ProcessLoadStoreTime;///< time of processing loads and stores
    static const char* ProcessCopyGepTime;	///< time of processing copys and geps
    static const char* UpdateCallGraphTime;	///< time of updating call graph

    static const char* TotalNumOfPointers;	///< Total PAG value node
    static const char* TotalNumOfObjects;	///< Total PAG object node
    static const char* TotalNumOfFieldObjects;	///< Total PAG field object node
    static const char* MaxStructSize;	///< Max struct size
    static const char* TotalNumOfEdges;	///< Total PAG edge number

    static const char* NumOfAddrs;		///< PAG addr edge
    static const char* NumOfLoads;		///< PAG load edge
    static const char* NumOfStores;		///< PAG store edge
    static const char* NumOfCopys;		///< PAG copy edge
    static const char* NumOfGeps;		///< PAG gep edge
    static const char* NumOfCalls;		///< PAG call edge
    static const char* NumOfReturns;	///< PAG return edge

    static const char* NumOfProcessedAddrs;		///< PAG processed addr edge
    static const char* NumOfProcessedLoads;		///< PAG processed load edge
    static const char* NumOfProcessedStores;	///< PAG processed store edge
    static const char* NumOfProcessedCopys;		///< PAG processed copy edge
    static const char* NumOfProcessedGeps;		///< PAG processed gep edge

    static const char* NumOfSfr;                ///< num of field representatives
    static const char* NumOfFieldExpand;

    static const char* NumOfFunctionObjs;	///< Function numbers
    static const char* NumOfGlobalObjs;	///< PAG global object node
    static const char* NumOfHeapObjs;	///< PAG heap object node
    static const char* NumOfStackObjs;	///< PAG stack object node

    static const char* NumOfObjsHasVarStruct;	///< PAG object node has var struct (maybe nested with array)
    static const char* NumOfObjsHasVarArray;	///< PAG object node has var array (maybe nested with struct)
    static const char* NumOfObjsHasConstStruct;	///< PAG object node has const struct (maybe nested with array)
    static const char* NumOfObjsHasConstArray;	///< PAG object node has const array (maybe nested with struct)
    static const char* NumOfNonPtrObjs;	///< PAG object node which is non pointer type object (do not have pts)
    static const char* NumOfConstantObjs;	///< PAG object node which is purely constant

    static const char* NumberOfFieldInsensitiveObj;
    static const char* NumberOfFieldSensitiveObj;

    static const char* NumOfPointers;	///< PAG value node, each of them maps to a llvm value
    static const char* NumOfGepFieldPointers;	///< PAG gep value node (field value, dynamically created dummy node)

    static const char* NumOfMemObjects;	///< PAG object node, each of them maps to a llvm value
    static const char* NumOfGepFieldObjects;	///< PAG gep object node (field obj, dynamically created dummy node)

    static const char* AveragePointsToSetSize;		///< Average points-to size of all variables
    static const char* AverageTopLevPointsToSetSize; ///< Average points-to size of top-level variables
    static const char* MaxPointsToSetSize;			///< Max points-to size

    static const char* NumOfIterations;	///< Number of iterations during resolution

    static const char* NumOfIndirectCallSites;	///< Number of indirect callsites
    static const char* NumOfIndirectEdgeSolved;	///< Number of indirect calledge resolved

    static const char* NumOfSCCDetection; ///< Number of scc detection performed
    static const char* NumOfCycles;   ///< Number of scc cycles detected
    static const char* NumOfPWCCycles;   ///< Number of scc cycles detected
    static const char* NumOfNodesInCycles; ///< Number of nodes in cycles detected
    static const char* MaxNumOfNodesInSCC;	///< max Number of nodes in one scc

    static const char* NumOfNullPointer;	///< Number of pointers points-to null

    typedef std::map<const char*,u32_t> NUMStatMap;

    typedef std::map<const char*,double> TIMEStatMap;

    PTAStat(PointerAnalysis* p);
    virtual ~PTAStat() {}

    virtual inline void startClk() {
        startTime = CLOCK_IN_MS();
    }
    virtual inline void endClk() {
        endTime = CLOCK_IN_MS();
    }
    static inline double getClk() {
        return CLOCK_IN_MS();
    }

    NUMStatMap generalNumMap;
    NUMStatMap PTNumStatMap;
    TIMEStatMap timeStatMap;
    NodeBS localVarInRecursion;

    double startTime;
    double endTime;

    virtual void performStat();

    virtual void printStat(string str = "");

    virtual void performStatPerQuery(NodeID ptr) {}

    virtual void printStatPerQuery(NodeID ptr, const PointsTo& pts) {}

    virtual void callgraphStat();
private:
    void bitcastInstStat();
    void branchStat();

    PointerAnalysis* pta;
    std::string moduleName;
};


/*
 * DDAStat.h
 *
 *  Created on: Sep 15, 2014
 *      Author: Yulei Sui
 */



class FlowDDA;
class ContextDDA;
class SVFG;
class PointerAnalysis;

/*!
 * Statistics of demand-driven analysis
 */
class DDAStat : public PTAStat {

public:
    DDAStat(FlowDDA* pta);
    DDAStat(ContextDDA* pta);

    u32_t _NumOfDPM;
    u32_t _NumOfStrongUpdates;
    u32_t _NumOfMustAliases;
    u32_t _NumOfInfeasiblePath;

    u64_t _NumOfStep;
    u64_t _NumOfStepInCycle;
    double _AnaTimePerQuery;
    double _AnaTimeCyclePerQuery;
    double _TotalTimeOfQueries;
    double _TotalTimeOfBKCondition;

    NodeBS _StrongUpdateStores;

    void performStatPerQuery(NodeID ptr);

    void performStat();

    void printStat();

    void printStatPerQuery(NodeID ptr, const PointsTo& pts);

    void getNumOfOOBQuery();

    inline void setMemUsageBefore(u32_t vmrss, u32_t vmsize) {
        _vmrssUsageBefore = vmrss;
        _vmsizeUsageBefore = vmsize;
    }

    inline void setMemUsageAfter(u32_t vmrss, u32_t vmsize) {
        _vmrssUsageAfter = vmrss;
        _vmsizeUsageAfter = vmsize;
    }

private:
    FlowDDA* flowDDA;
    ContextDDA* contextDDA;

    u32_t _TotalNumOfQuery;
    u32_t _TotalNumOfOutOfBudgetQuery;
    u32_t _TotalNumOfDPM;
    u32_t _TotalNumOfStrongUpdates;
    u32_t _TotalNumOfMustAliases;
    u32_t _TotalNumOfInfeasiblePath;

    u32_t _TotalNumOfStep;
    u32_t _TotalNumOfStepInCycle;

    u32_t _NumOfIndCallEdgeSolved;
    u32_t _MaxCPtsSize;
    u32_t _MaxPtsSize;
    u32_t _TotalCPtsSize;
    u32_t _TotalPtsSize;
    u32_t _NumOfNullPtr;
    u32_t _NumOfConstantPtr;
    u32_t _NumOfBlackholePtr;

    u32_t _vmrssUsageBefore;
    u32_t _vmrssUsageAfter;
    u32_t _vmsizeUsageBefore;
    u32_t _vmsizeUsageAfter;

    double _AvgNumOfDPMAtSVFGNode;
    u32_t _MaxNumOfDPMAtSVFGNode;

    NUMStatMap NumPerQueryStatMap;

    void initDefault();

public:
    SVFG* getSVFG() const;

    PointerAnalysis* getPTA() const;

    inline NodeBS& getStrongUpdateStores() {
        return _StrongUpdateStores;
    }
};

/*!
 * Pointer analysis implementation which uses conditional points-to map data structure (context/path sensitive analysis)
 */
template<class Cond>
class CondPTAImpl : public PointerAnalysis {

public:
    typedef CondVar<Cond> CVar;
    typedef CondStdSet<CVar>  CPtSet;
    typedef PTData<CVar,CPtSet> PTDataTy;	         /// Points-to data structure type
    typedef std::map<NodeID,PointsTo> PtrToBVPtsMap; /// map a pointer to its BitVector points-to representation
    typedef std::map<NodeID,CPtSet> PtrToCPtsMap;	 /// map a pointer to its conditional points-to set

    /// Constructor
    CondPTAImpl(PointerAnalysis::PTATY type) : PointerAnalysis(type), normalized(false) {
        if (type == PathS_DDA || type == Cxt_DDA)
            ptD = new PTDataTy();
        else
            assert(false && "no points-to data available");
    }

    /// Destructor
    virtual ~CondPTAImpl() {
        destroy();
    }

    /// Release memory
    inline void destroy() {
        delete ptD;
        ptD = NULL;
    }

    /// Get points-to data
    inline PTDataTy* getPTDataTy() const {
        return ptD;
    }

    /// Get points-to and reverse points-to
    ///@{
    virtual inline CPtSet& getPts(CVar id) {
        return ptD->getPts(id);
    }
    virtual inline CPtSet& getRevPts(CVar nodeId) {
        return ptD->getRevPts(nodeId);
    }
    //@}

    /// Clear all data
    virtual inline void clearPts() {
        ptD->clear();
    }

    /// Whether cpts1 and cpts2 have overlap points-to targets
    bool overlap(const CPtSet& cpts1, const CPtSet& cpts2) const {
        for (typename CPtSet::const_iterator it1 = cpts1.begin(); it1 != cpts1.end(); ++it1) {
            for (typename CPtSet::const_iterator it2 = cpts2.begin(); it2 != cpts2.end(); ++it2) {
                if(isSameVar(*it1,*it2))
                    return true;
            }
        }
        return false;
    }

    /// Expand all fields of an aggregate in all points-to sets
    void expandFIObjs(const CPtSet& cpts, CPtSet& expandedCpts) {
        expandedCpts = cpts;;
        for(typename CPtSet::const_iterator cit = cpts.begin(), ecit=cpts.end(); cit!=ecit; ++cit) {
            if(pag->getBaseObjNode(cit->get_id())==cit->get_id()) {
                NodeBS& fields = pag->getAllFieldsObjNode(cit->get_id());
                for(NodeBS::iterator it = fields.begin(), eit = fields.end(); it!=eit; ++it) {
                    CVar cvar(cit->get_cond(),*it);
                    expandedCpts.set(cvar);
                }
            }
        }
    }

protected:

    /// Finalization of pointer analysis, and normalize points-to information to Bit Vector representation
    virtual void finalize() {
        NormalizePointsTo();
        PointerAnalysis::finalize();
    }
    /// Union/add points-to, and add the reverse points-to for node collapse purpose
    /// To be noted that adding reverse pts might incur 10% total overhead during solving
    //@{
    virtual inline bool unionPts(CVar id, const CPtSet& target) {
        return ptD->unionPts(id, target);
    }

    virtual inline bool unionPts(CVar id, CVar ptd) {
        return ptD->unionPts(id,ptd);
    }

    virtual inline bool addPts(CVar id, CVar ptd) {
        return ptD->addPts(id,ptd);
    }
    //@}

    /// Internal interface to be used for conditional points-to set queries
    //@{
    inline bool mustAlias(const CVar& var1, const CVar& var2) {
        if(isSameVar(var1,var2))
            return true;

        bool singleton = !(isHeapMemObj(var1.get_id()) || isLocalVarInRecursiveFun(var1.get_id()));
        if(isCondCompatible(var1.get_cond(),var2.get_cond(),singleton) == false)
            return false;

        const CPtSet& cpts1 = getPts(var1);
        const CPtSet& cpts2 = getPts(var2);
        return (contains(cpts1,cpts2) && contains(cpts2,cpts1));
    }

    //  Whether cpts1 contains all points-to targets of pts2
    bool contains(const CPtSet& cpts1, const CPtSet& cpts2) {
        if (cpts1.empty() || cpts2.empty())
            return false;

        for (typename CPtSet::const_iterator it2 = cpts2.begin(); it2 != cpts2.end(); ++it2) {
            bool hasObj = false;
            for (typename CPtSet::const_iterator it1 = cpts1.begin(); it1 != cpts1.end(); ++it1) {
                if(isSameVar(*it1,*it2)) {
                    hasObj = true;
                    break;
                }
            }
            if(hasObj == false)
                return false;
        }
        return true;
    }

    /// Whether two pointers/objects are the same one by considering their conditions
    bool isSameVar(const CVar& var1, const CVar& var2) const {
        if(var1.get_id() != var2.get_id())
            return false;

        /// we distinguish context sensitive memory allocation here
        bool singleton = !(isHeapMemObj(var1.get_id()) || isLocalVarInRecursiveFun(var1.get_id()));
        return isCondCompatible(var1.get_cond(),var2.get_cond(),singleton);
    }
    //@}

    /// Normalize points-to information to BitVector/conditional representation
    virtual void NormalizePointsTo() {
        normalized = true;
        const typename PTDataTy::PtsMap& ptsMap = getPTDataTy()->getPtsMap();
        for(typename PTDataTy::PtsMap::const_iterator it = ptsMap.begin(), eit=ptsMap.end(); it!=eit; ++it) {
            for(typename CPtSet::const_iterator cit = it->second.begin(), ecit=it->second.end(); cit!=ecit; ++cit) {
                ptrToBVPtsMap[(it->first).get_id()].set(cit->get_id());
                objToBVRevPtsMap[cit->get_id()].set((it->first).get_id());
                ptrToCPtsMap[(it->first).get_id()].set(*cit);
            }
        }
    }
    /// Points-to data
    PTDataTy* ptD;
    /// Normalized flag
    bool normalized;
    /// Normal points-to representation (without conditions)
    PtrToBVPtsMap ptrToBVPtsMap;
    /// Normal points-to representation (without conditions)
    PtrToBVPtsMap objToBVRevPtsMap;
    /// Conditional points-to representation (with conditions)
    PtrToCPtsMap ptrToCPtsMap;
public:
    /// Print out conditional pts
    virtual void dumpCPts() {
        ptD->dumpPTData();
    }
    /// Given a conditional pts return its bit vector points-to
    virtual inline PointsTo getBVPointsTo(const CPtSet& cpts) const {
        PointsTo pts;
        for(typename CPtSet::const_iterator cit = cpts.begin(), ecit=cpts.end(); cit!=ecit; ++cit)
            pts.set(cit->get_id());
        return pts;
    }
    /// Given a pointer return its bit vector points-to
    virtual inline PointsTo& getPts(NodeID ptr) {
        assert(normalized && "Pts of all context-var have to be merged/normalized. Want to use getPts(CVar cvar)??");
        return ptrToBVPtsMap[ptr];
    }
    /// Given a pointer return its conditional points-to
    virtual inline const CPtSet& getCondPointsTo(NodeID ptr) {
        assert(normalized && "Pts of all context-vars have to be merged/normalized. Want to use getPts(CVar cvar)??");
        return ptrToCPtsMap[ptr];
    }
    /// Given an object return all pointers points to this object
    virtual inline PointsTo& getRevPts(NodeID obj) {
        assert(normalized && "Pts of all context-var have to be merged/normalized. Want to use getPts(CVar cvar)??");
        return objToBVRevPtsMap[obj];
    }

    /// Interface expose to users of our pointer analysis, given Location infos
    virtual inline AliasResult alias(const MemoryLocation &LocA,
                                           const MemoryLocation  &LocB) {
        return alias(LocA.Ptr, LocB.Ptr);
    }
    /// Interface expose to users of our pointer analysis, given Value infos
    virtual inline AliasResult alias(const Value* V1, const Value* V2) {
        return  alias(pag->getValueNode(V1),pag->getValueNode(V2));
    }
    /// Interface expose to users of our pointer analysis, given two pointers
    virtual inline AliasResult alias(NodeID node1, NodeID node2) {
        return alias(getCondPointsTo(node1),getCondPointsTo(node2));
    }
    /// Interface expose to users of our pointer analysis, given conditional variables
    virtual AliasResult alias(const CVar& var1, const CVar& var2) {
        return alias(getPts(var1),getPts(var2));
    }
    /// Interface expose to users of our pointer analysis, given two conditional points-to sets
    virtual inline AliasResult alias(const CPtSet& pts1, const CPtSet& pts2) {
        CPtSet cpts1;
        expandFIObjs(pts1,cpts1);
        CPtSet cpts2;
        expandFIObjs(pts2,cpts2);
        if (containBlackHoleNode(cpts1) || containBlackHoleNode(cpts2))
            return llvm::MayAlias;
        else if(this->getAnalysisTy()==PathS_DDA && contains(cpts1,cpts2) && contains(cpts2,cpts1)) {
            return llvm::MustAlias;
        }
        else if(overlap(cpts1,cpts2))
            return llvm::MayAlias;
        else
            return llvm::NoAlias;
    }
    /// Test blk node for cpts
    inline bool containBlackHoleNode(const CPtSet& cpts) {
        for(typename CPtSet::const_iterator cit = cpts.begin(), ecit=cpts.end(); cit!=ecit; ++cit) {
            if(cit->get_id() == pag->getBlackHoleNode())
                return true;
        }
        return false;
    }
    /// Test constant node for cpts
    inline bool containConstantNode(const CPtSet& cpts) {
        for(typename CPtSet::const_iterator cit = cpts.begin(), ecit=cpts.end(); cit!=ecit; ++cit) {
            if(cit->get_id() == pag->getConstantNode())
                return true;
        }
        return false;
    }
    /// Whether two conditions are compatible (to be implemented by child class)
    virtual bool isCondCompatible(const Cond& cxt1, const Cond& cxt2, bool singleton) const = 0;

    /// Dump points-to information of top-level pointers
    void dumpTopLevelPtsTo() {
        for (NodeSet::iterator nIter = this->getAllValidPtrs().begin(); nIter != this->getAllValidPtrs().end(); ++nIter) {
            const PAGNode* node = this->getPAG()->getPAGNode(*nIter);
            if (this->getPAG()->isValidTopLevelPtr(node)) {
                if (SVFUtil::isa<DummyObjPN>(node)) {
                    SVFUtil::outs() << "##<Blackhole or constant> id:" << node->getId();
                }
                else if (!SVFUtil::isa<DummyValPN>(node)) {
                    SVFUtil::outs() << "##<" << node->getValue()->getName() << "> ";
                    SVFUtil::outs() << "Source Loc: " << SVFUtil::getSourceLoc(node->getValue());
                }

                const PointsTo& pts = getPts(node->getId());
                SVFUtil::outs() << "\nNodeID " << node->getId() << " ";
                if (pts.empty()) {
                    SVFUtil::outs() << "\t\tPointsTo: {empty}\n\n";
                } else {
                    SVFUtil::outs() << "\t\tPointsTo: { ";
                    for (PointsTo::iterator it = pts.begin(), eit = pts.end(); it != eit; ++it)
                        SVFUtil::outs() << *it << " ";
                    SVFUtil::outs() << "}\n\n";
                }
            }
        }
    }
};
/*
 * ContextDDA.h
 *
 *  Created on: Aug 17, 2014
 *      Author: Yulei Sui
 */



class FlowDDA;
class DDAClient;

/*!
 * Context-, Flow- Sensitive Demand-driven Analysis
 */

%template(CondPTAImplContext) CondPTAImpl<ContextCond>;
%template(DDAVFSolverCCC) DDAVFSolver<CxtVar,CxtPtSet,CxtLocDPItem>;
class ContextDDA : public CondPTAImpl<ContextCond>, public DDAVFSolver<CxtVar,CxtPtSet,CxtLocDPItem> {

public:
    /// Constructor
    ContextDDA(SVFModule mod, DDAClient* client);

    /// Destructor
    virtual ~ContextDDA();

    /// Initialization of the analysis
    virtual void initialize(SVFModule module);

    /// Finalize analysis
    virtual inline void finalize() {
        CondPTAImpl<ContextCond>::finalize();
    }

    /// dummy analyze method
    virtual void analyze(SVFModule mod) {}

    /// Compute points-to set for an unconditional pointer
    void computeDDAPts(NodeID id);

    /// Compute points-to set for a context-sensitive pointer
    const CxtPtSet& computeDDAPts(const CxtVar& cxtVar);

    /// Handle out-of-budget dpm
    void handleOutOfBudgetDpm(const CxtLocDPItem& dpm);

    /// Override parent method
    CxtPtSet getConservativeCPts(const CxtLocDPItem& dpm) {
        const PointsTo& pts =  getAndersenAnalysis()->getPts(dpm.getCurNodeID());
        CxtPtSet tmpCPts;
        ContextCond cxt;
        for (PointsTo::iterator piter = pts.begin(); piter != pts.end(); ++piter) {
            CxtVar var(cxt,*piter);
            tmpCPts.set(var);
        }
        return tmpCPts;
    }

    /// Override parent method
    virtual inline NodeID getPtrNodeID(const CxtVar& var) const {
        return var.get_id();
    }
    /// Handle condition for context or path analysis (backward analysis)
    virtual bool handleBKCondition(CxtLocDPItem& dpm, const SVFGEdge* edge);

    /// we exclude concrete heap given the following conditions:
    /// (1) concrete calling context (not involved in recursion and not exceed the maximum context limit)
    /// (2) not inside loop
    bool isHeapCondMemObj(const CxtVar& var, const StoreSVFGNode* store);

    /// refine indirect call edge
    bool testIndCallReachability(CxtLocDPItem& dpm, const Function* callee, CallSite cs);

    /// get callsite id from call, return 0 if it is a spurious call edge
    CallSiteID getCSIDAtCall(CxtLocDPItem& dpm, const SVFGEdge* edge);

    /// get callsite id from return, return 0 if it is a spurious return edge
    CallSiteID getCSIDAtRet(CxtLocDPItem& dpm, const SVFGEdge* edge);


    /// Pop recursive callsites
    inline virtual void popRecursiveCallSites(CxtLocDPItem& dpm) {
        ContextCond& cxtCond = dpm.getCond();
        cxtCond.setNonConcreteCxt();
        CallStrCxt& cxt = cxtCond.getContexts();
        while(!cxt.empty() && isEdgeInRecursion(cxt.back())) {
            cxt.pop_back();
        }
    }
    /// Whether call/return inside recursion
    inline virtual bool isEdgeInRecursion(CallSiteID csId) {
        const Function* caller = getPTACallGraph()->getCallerOfCallSite(csId);
        const Function* callee = getPTACallGraph()->getCalleeOfCallSite(csId);
        return inSameCallGraphSCC(caller, callee);
    }
    /// Update call graph.
    //@{
    void updateCallGraphAndSVFG(const CxtLocDPItem& dpm,CallSite cs,SVFGEdgeSet& svfgEdges)
    {
        CallEdgeMap newEdges;
        resolveIndCalls(cs, getBVPointsTo(getCachedPointsTo(dpm)), newEdges);
        for (CallEdgeMap::const_iterator iter = newEdges.begin(),eiter = newEdges.end(); iter != eiter; iter++) {
            CallSite newcs = iter->first;
            const FunctionSet & functions = iter->second;
            for (FunctionSet::const_iterator func_iter = functions.begin(); func_iter != functions.end(); func_iter++) {
                const Function * func = *func_iter;
                getSVFG()->connectCallerAndCallee(newcs, func, svfgEdges);
            }
        }
    }
    //@}

    /// Return TRUE if this edge is inside a SVFG SCC, i.e., src node and dst node are in the same SCC on the SVFG.
    inline bool edgeInCallGraphSCC(const SVFGEdge* edge) {
        const BasicBlock* srcBB = edge->getSrcNode()->getBB();
        const BasicBlock* dstBB = edge->getDstNode()->getBB();

        if(srcBB && dstBB)
            return inSameCallGraphSCC(srcBB->getParent(),dstBB->getParent());

        assert(edge->isRetVFGEdge() == false && "should not be an inter-procedural return edge" );

        return false;
    }

    /// processGep node
    CxtPtSet processGepPts(const GepSVFGNode* gep, const CxtPtSet& srcPts);

    /// Handle Address SVFGNode to add proper conditional points-to
    void handleAddr(CxtPtSet& pts,const CxtLocDPItem& dpm,const AddrSVFGNode* addr) {
        NodeID srcID = addr->getPAGSrcNodeID();
        /// whether this object is set field-insensitive during pre-analysis
        if (isFieldInsensitive(srcID))
            srcID = getFIObjNode(srcID);

        CxtVar var(dpm.getCond(),srcID);
        addDDAPts(pts,var);
        DBOUT(DDDA, SVFUtil::outs() << "\t add points-to target " << var << " to dpm ");
        DBOUT(DDDA, dpm.dump());
    }

    /// Propagate along indirect value-flow if two objects of load and store are same
    virtual inline bool propagateViaObj(const CxtVar& storeObj, const CxtVar& loadObj) {
        return isSameVar(storeObj,loadObj);
    }

    /// Whether two call string contexts are compatible which may represent the same memory object
    /// compare with call strings from last few callsite ids (most recent ids to objects):
    /// compatible : (e.g., 123 == 123, 123 == 23). not compatible (e.g., 123 != 423)
    inline bool isCondCompatible(const ContextCond& cxt1, const ContextCond& cxt2, bool singleton) const;

    /// Whether this edge is treated context-insensitively
    bool isInsensitiveCallRet(const SVFGEdge* edge) {
        return insensitveEdges.find(edge) != insensitveEdges.end();
    }
    /// Return insensitive edge set
    inline ConstSVFGEdgeSet& getInsensitiveEdgeSet() {
        return insensitveEdges;
    }
    /// dump context call strings
    virtual inline void dumpContexts(const ContextCond& cxts) {
        SVFUtil::outs() << cxts.toString() << "\n";
    }

    virtual const std::string PTAName() const {
        return "Context Sensitive DDA";
    }

private:
    ConstSVFGEdgeSet insensitveEdges;///< insensitive call-return edges
    FlowDDA* flowDDA;			///< downgrade to flowDDA if out-of-budget
    DDAClient* _client;			///< DDA client
    PTACFInfoBuilder loopInfoBuilder; ///< LoopInfo
};


/*
 * @file: DDAClient.h
 * @author: yesen
 * @date: 4 Feb 2015
 *
 * LICENSE
 *
 */

/**
 * General DDAClient which queries all top level pointers by default.
 */
class DDAClient {
public:
    DDAClient(SVFModule mod) : pag(NULL), module(mod), curPtr(0), solveAll(true) {}

    virtual ~DDAClient() {}

    virtual inline void initialise(SVFModule module) {}

    /// Collect candidate pointers for query.
    virtual inline NodeSet& collectCandidateQueries(PAG* p) {
        setPAG(p);
        if (solveAll)
            candidateQueries = pag->getAllValidPtrs();
        else {
            for (NodeSet::iterator it = userInput.begin(), eit = userInput.end(); it != eit; ++it)
                addCandidate(*it);
        }
        return candidateQueries;
    }
    /// Get candidate queries
    inline const NodeSet& getCandidateQueries() const {
        return candidateQueries;
    }

    /// Call back used by DDAVFSolver.
    virtual inline void handleStatement(const SVFGNode* stmt, NodeID var) {}
    /// Set PAG graph.
    inline void setPAG(PAG* g) {
        pag = g;
    }
    /// Set the pointer being queried.
    void setCurrentQueryPtr(NodeID ptr) {
        curPtr = ptr;
    }
    /// Set pointer to be queried by DDA analysis.
    void setQuery(NodeID ptr) {
        userInput.insert(ptr);
        solveAll = false;
    }
    /// Get LLVM module
    inline SVFModule getModule() const {
        return module;
    }
    virtual void answerQueries(PointerAnalysis* pta);

    virtual inline void performStat(PointerAnalysis* pta) {}

    virtual inline void collectWPANum(SVFModule mod) {}
protected:
    void addCandidate(NodeID id) {
        if (pag->isValidTopLevelPtr(pag->getPAGNode(id)))
            candidateQueries.insert(id);
    }

    PAG*   pag;					///< PAG graph used by current DDA analysis
    SVFModule module;		///< LLVM module
    NodeID curPtr;				///< current pointer being queried
    NodeSet candidateQueries;	///< store all candidate pointers to be queried

private:
    NodeSet userInput;           ///< User input queries
    bool solveAll;				///< TRUE if all top level pointers are being queried
};


/**
 * DDA client with function pointers as query candidates.
 */
class FunptrDDAClient : public DDAClient {
private:
    typedef std::map<NodeID,CallSite> VTablePtrToCallSiteMap;
    VTablePtrToCallSiteMap vtableToCallSiteMap;
public:
    FunptrDDAClient(SVFModule module) : DDAClient(module) {}
    ~FunptrDDAClient() {}

    /// Only collect function pointers as query candidates.
    virtual inline NodeSet& collectCandidateQueries(PAG* p) {
        setPAG(p);
        for(PAG::CallSiteToFunPtrMap::const_iterator it = pag->getIndirectCallsites().begin(),
                eit = pag->getIndirectCallsites().end(); it!=eit; ++it) {
            if (cppUtil::isVirtualCallSite(it->first)) {
                const Value *vtblPtr = cppUtil::getVCallVtblPtr(it->first);
                assert(pag->hasValueNode(vtblPtr) && "not a vtable pointer?");
                NodeID vtblId = pag->getValueNode(vtblPtr);
                addCandidate(vtblId);
                vtableToCallSiteMap[vtblId] = it->first;
            } else {
                addCandidate(it->second);
            }
        }
        return candidateQueries;
    }
    virtual void performStat(PointerAnalysis* pta);
};



//===- BasicTypes.h -- Basic types used in SVF-------------------------------//
//
//                     SVF: Static Value-Flow Analysis
//
// Copyright (C) <2013-2017>  <Yulei Sui>
//

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//===----------------------------------------------------------------------===//

/*
 * BasicTypes.h
 *
 *  Created on: Apr 1, 2014
 *      Author: Yulei Sui
 */








/// LLVM debug macros, define type of your DEBUG model of each pass
#define DBOUT(TYPE, X) 	DEBUG_WITH_TYPE(TYPE, X)
#define DOSTAT(X) 	X
#define DOTIMESTAT(X) 	X

/// General debug flag is for each phase of a pass, it is often in a colorful output format
#define DGENERAL "general"

#define DPAGBuild "pag"
#define DMemModel "mm"
#define DMemModelCE "mmce"
#define DCOMModel "comm"
#define DDDA "dda"
#define DDumpPT "dumppt"
#define DRefinePT "sbpt"
#define DCache "cache"
#define DWPA "wpa"
#define DMSSA "mssa"
#define DInstrument "ins"
#define DAndersen "ander"
#define DSaber "saber"
#define DMTA "mta"
#define DCHA "cha"

/*
 * Number of clock ticks per second. A clock tick is the unit by which
 * processor time is measured and is returned by 'clock'.
 */
#define TIMEINTERVAL 1000
#define CLOCK_IN_MS() (clock() / (CLOCKS_PER_SEC / TIMEINTERVAL))

class BddCond;


/// LLVM Basic classes

/// LLVM outputs

/// LLVM types

/// LLVM data layout

/// LLVM Aliases and constants

/// LLVM metadata


/// LLVM Instructions

/// LLVM scalar evolution

/// LLVM Dominators

/// LLVM Iterators
//===- PTAStat.h -- Base class for statistics---------------------------------//
//
//                     SVF: Static Value-Flow Analysis
//
// Copyright (C) <2013-2017>  <Yulei Sui>
//

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//===----------------------------------------------------------------------===//

/*
 * AndersenStat.h
 *
 *  Created on: Oct 12, 2013
 *      Author: Yulei Sui
 */



using namespace std;



//===- SVFUtil.h -- Analysis helper functions----------------------------//
//
//                     SVF: Static Value-Flow Analysis
//
// Copyright (C) <2013-2017>  <Yulei Sui>
//

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//===----------------------------------------------------------------------===//

/*
 * SVFUtil.h
 *
 *  Created on: Apr 11, 2013
 *      Author: Yulei Sui, dye
 */



/*
 * Util class to assist pointer analysis
 */
namespace SVFUtil {

/// Overwrite llvm::outs()
inline raw_ostream &outs(){
	return llvm::outs();
}

/// Overwrite llvm::errs()
inline raw_ostream &errs(){
	return llvm::errs();
}

/// Return true if this function is llvm dbg intrinsic function/instruction
//@{
inline bool isIntrinsicDbgFun(const Function* fun) {
    return fun->getName().startswith("llvm.dbg.declare") ||
           fun->getName().startswith("llvm.dbg.value");
}
bool isInstrinsicDbgInst(const Instruction* inst);
//@}

/// This function servers a allocation wrapper detector
inline bool isAnAllocationWraper(const Instruction *inst) {
    return false;
}
/// Whether an instruction is a call or invoke instruction
inline bool isCallSite(const Instruction* inst) {
    return SVFUtil::isa<CallInst>(inst) || SVFUtil::isa<InvokeInst>(inst);
}
/// Whether an instruction is a callsite in the application code, excluding llvm intrinsic calls
inline bool isNonInstricCallSite(const Instruction* inst) {
	if(isInstrinsicDbgInst(inst))
		return false;
    return isCallSite(inst);
}
/// Whether an instruction is a return instruction
inline bool isReturn(const Instruction* inst) {
    return SVFUtil::isa<ReturnInst>(inst);
}
/// Return LLVM function if this value is
inline const Function* getLLVMFunction(const Value* val) {
    const Function *fun = SVFUtil::dyn_cast<Function>(val->stripPointerCasts());
    return fun;
}
/// Return LLVM callsite given a instruction
inline CallSite getLLVMCallSite(const Instruction* inst) {
    assert(SVFUtil::isa<CallInst>(inst)|| SVFUtil::isa<InvokeInst>(inst));
    CallSite cs(const_cast<Instruction*>(inst));
    return cs;
}

/// Get the definition of a function across multiple modules
inline const Function* getDefFunForMultipleModule(const Function* fun) {
	if(fun == NULL) return NULL;

    SVFModule svfModule;
    if (fun->isDeclaration() && svfModule.hasDefinition(fun))
        fun = svfModule.getDefinition(fun);
    return fun;
}

/// Return callee of a callsite. Return null if this is an indirect call
//@{
inline const Function* getCallee(const CallSite cs) {
    // FIXME: do we need to strip-off the casts here to discover more library functions
    Function *callee = SVFUtil::dyn_cast<Function>(cs.getCalledValue()->stripPointerCasts());
    return getDefFunForMultipleModule(callee);
}

inline const Function* getCallee(const Instruction *inst) {
    if (!SVFUtil::isa<CallInst>(inst) && !SVFUtil::isa<InvokeInst>(inst))
        return NULL;
    CallSite cs(const_cast<Instruction*>(inst));
    return getCallee(cs);
}
//@}

/// Return true if the call is an external call (external library in function summary table)
//@{
inline bool isExtCall(const Function* fun) {
    return fun && ExtAPI::getExtAPI()->is_ext(fun);
}

inline bool isExtCall(const CallSite cs) {
    return isExtCall(getCallee(cs));
}

inline bool isExtCall(const Instruction *inst) {
    return isExtCall(getCallee(inst));
}
//@}

/// Return true if the call is a heap allocator/reallocator
//@{
/// note that these two functions are not suppose to be used externally
inline bool isHeapAllocExtFunViaRet(const Function *fun) {
    return fun && (ExtAPI::getExtAPI()->is_alloc(fun)
                   || ExtAPI::getExtAPI()->is_realloc(fun));
}
inline bool isHeapAllocExtFunViaArg(const Function *fun) {
    return fun && ExtAPI::getExtAPI()->is_arg_alloc(fun);
}

/// interfaces to be used externally
inline bool isHeapAllocExtCallViaRet(const CallSite cs) {
    bool isPtrTy = cs.getInstruction()->getType()->isPointerTy();
    return isPtrTy && isHeapAllocExtFunViaRet(getCallee(cs));
}

inline bool isHeapAllocExtCallViaRet(const Instruction *inst) {
    bool isPtrTy = inst->getType()->isPointerTy();
    return isPtrTy && isHeapAllocExtFunViaRet(getCallee(inst));
}

inline bool isHeapAllocExtCallViaArg(const CallSite cs) {
    return isHeapAllocExtFunViaArg(getCallee(cs));
}

inline bool isHeapAllocExtCallViaArg(const Instruction *inst) {
    return isHeapAllocExtFunViaArg(getCallee(inst));
}

inline bool isHeapAllocExtCall(const CallSite cs) {
    return isHeapAllocExtCallViaRet(cs) || isHeapAllocExtCallViaArg(cs);
}

inline bool isHeapAllocExtCall(const Instruction *inst) {
    return isHeapAllocExtCallViaRet(inst) || isHeapAllocExtCallViaArg(inst);
}
//@}

/// Get the position of argument that holds an allocated heap object.
//@{
inline int getHeapAllocHoldingArgPosition(const Function *fun) {
    return ExtAPI::getExtAPI()->get_alloc_arg_pos(fun);
}

inline int getHeapAllocHoldingArgPosition(const CallSite cs) {
    return getHeapAllocHoldingArgPosition(getCallee(cs));
}

inline int getHeapAllocHoldingArgPosition(const Instruction *inst) {
    return getHeapAllocHoldingArgPosition(getCallee(inst));
}
//@}

/// Return true if the call is a heap reallocator
//@{
/// note that this function is not suppose to be used externally
inline bool isReallocExtFun(const Function *fun) {
    return fun && (ExtAPI::getExtAPI()->is_realloc(fun));
}

inline bool isReallocExtCall(const CallSite cs) {
    bool isPtrTy = cs.getInstruction()->getType()->isPointerTy();
    return isPtrTy && isReallocExtFun(getCallee(cs));
}

inline bool isReallocExtCall(const Instruction *inst) {
    bool isPtrTy = inst->getType()->isPointerTy();
    return isPtrTy && isReallocExtFun(getCallee(inst));
}
//@}

/// Return true if the call is a heap dealloc or not
//@{
/// note that this function is not suppose to be used externally
inline bool isDeallocExtFun(const Function *fun) {
    return fun && (ExtAPI::getExtAPI()->is_dealloc(fun));
}

inline bool isDeallocExtCall(const CallSite cs) {
    return isDeallocExtFun(getCallee(cs));
}

inline bool isDeallocExtCall(const Instruction *inst) {
    return isDeallocExtFun(getCallee(inst));
}
//@}


/// Return true if the call is a static global call
//@{
/// note that this function is not suppose to be used externally
inline bool isStaticExtFun(const Function *fun) {
    return fun && ExtAPI::getExtAPI()->has_static(fun);
}

inline bool isStaticExtCall(const CallSite cs) {
    bool isPtrTy = cs.getInstruction()->getType()->isPointerTy();
    return isPtrTy && isStaticExtFun(getCallee(cs));
}

inline bool isStaticExtCall(const Instruction *inst) {
    bool isPtrTy = inst->getType()->isPointerTy();
    return isPtrTy && isStaticExtFun(getCallee(inst));
}
//@}

/// Return true if the call is a static global call
//@{
inline bool isHeapAllocOrStaticExtCall(const CallSite cs) {
    return isStaticExtCall(cs) || isHeapAllocExtCall(cs);
}

inline bool isHeapAllocOrStaticExtCall(const Instruction *inst) {
    return isStaticExtCall(inst) || isHeapAllocExtCall(inst);
}
//@}

/// Return external call type
inline ExtAPI::extf_t extCallTy(const Function* fun) {
    return ExtAPI::getExtAPI()->get_type(fun);
}

/// Get the reference type of heap/static object from an allocation site.
//@{
inline const PointerType *getRefTypeOfHeapAllocOrStatic(const CallSite cs) {
    const PointerType *refType = NULL;
    // Case 1: heap object held by *argument, we should get its element type.
    if (isHeapAllocExtCallViaArg(cs)) {
        int argPos = getHeapAllocHoldingArgPosition(cs);
        const Value *arg = cs.getArgument(argPos);
        if (const PointerType *argType = SVFUtil::dyn_cast<PointerType>(arg->getType()))
            refType = SVFUtil::dyn_cast<PointerType>(argType->getElementType());
    }
    // Case 2: heap/static object held by return value.
    else {
        assert((isStaticExtCall(cs) || isHeapAllocExtCallViaRet(cs))
               && "Must be heap alloc via ret, or static allocation site");
        refType = SVFUtil::dyn_cast<PointerType>(cs.getType());
    }
    assert(refType && "Allocated object must be held by a pointer-typed value.");
    return refType;
}

inline const PointerType *getRefTypeOfHeapAllocOrStatic(const Instruction *inst) {
    CallSite cs(const_cast<Instruction*>(inst));
    return getRefTypeOfHeapAllocOrStatic(cs);
}
//@}

/// Return true if this is a thread creation call
///@{
inline bool isThreadForkCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isTDFork(cs);
}
inline bool isThreadForkCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isTDFork(inst);
}
//@}

/// Return true if this is a hare_parallel_for call
///@{
inline bool isHareParForCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isHareParFor(cs);
}
inline bool isHareParForCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isHareParFor(inst);
}
//@}

/// Return true if this is a thread join call
///@{
inline bool isThreadJoinCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isTDJoin(cs);
}
inline bool isThreadJoinCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isTDJoin(inst);
}
//@}

/// Return true if this is a thread exit call
///@{
inline bool isThreadExitCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isTDExit(cs);
}
inline bool isThreadExitCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isTDExit(inst);
}
//@}

/// Return true if this is a lock acquire call
///@{
inline bool isLockAquireCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isTDAcquire(cs);
}
inline bool isLockAquireCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isTDAcquire(inst);
}
//@}

/// Return true if this is a lock acquire call
///@{
inline bool isLockReleaseCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isTDRelease(cs);
}
inline bool isLockReleaseCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isTDRelease(inst);
}
//@}

/// Return true if this is a barrier wait call
//@{
inline bool isBarrierWaitCall(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->isTDBarWait(cs);
}
inline bool isBarrierWaitCall(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->isTDBarWait(inst);
}
//@}

/// Return thread fork function
//@{
inline const Value* getForkedFun(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->getForkedFun(cs);
}
inline const Value* getForkedFun(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->getForkedFun(inst);
}
//@}

/// Return sole argument of the thread routine
//@{
inline const Value* getActualParmAtForkSite(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->getActualParmAtForkSite(cs);
}
inline const Value* getActualParmAtForkSite(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->getActualParmAtForkSite(inst);
}
//@}

/// Return the task function of the parallel_for routine
//@{
inline const Value* getTaskFuncAtHareParForSite(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->getTaskFuncAtHareParForSite(cs);
}
inline const Value* getTaskFuncAtHareParForSite(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->getTaskFuncAtHareParForSite(inst);
}
//@}

/// Return the task data argument of the parallel_for rountine
//@{
inline const Value* getTaskDataAtHareParForSite(const CallSite cs) {
    return ThreadAPI::getThreadAPI()->getTaskDataAtHareParForSite(cs);
}
inline const Value* getTaskDataAtHareParForSite(const Instruction *inst) {
    return ThreadAPI::getThreadAPI()->getTaskDataAtHareParForSite(inst);
}
//@}

/// Return true if this value refers to a object
bool isObject (const Value * ref);


/// Method for dead function, which does not have any possible caller
/// function address is not taken and never be used in call or invoke instruction
//@{
/// whether this is a function without any possible caller?
bool isDeadFunction (const Function * fun);

/// whether this is an argument in dead function
inline bool ArgInDeadFunction (const Value * val) {
    return SVFUtil::isa<Argument>(val)
           && isDeadFunction(SVFUtil::cast<Argument>(val)->getParent());
}
//@}

/// Program entry function e.g. main
//@{
/// Return true if this is a program entry function (e.g. main)
inline bool isProgEntryFunction (const Function * fun) {
    return fun && fun->getName().str() == "main";
}

/// Get program entry function from module.
inline const Function* getProgEntryFunction(SVFModule svfModule) {
    for (SVFModule::const_iterator it = svfModule.begin(), eit = svfModule.end(); it != eit; ++it) {
        const Function *fun = *it;
        if (isProgEntryFunction(fun))
            return (fun);
    }
    return NULL;
}

/// Return true if this is an argument of a program entry function (e.g. main)
inline bool ArgInProgEntryFunction (const Value * val) {
    return SVFUtil::isa<Argument>(val)
           && isProgEntryFunction(SVFUtil::cast<Argument>(val)->getParent());
}
/// Return true if this is value in a dead function (function without any caller)
bool isPtrInDeadFunction (const Value * value);
//@}

/// Return true if this is a program exit function call
//@{
inline bool isProgExitFunction (const Function * fun) {
    return fun && (fun->getName().str() == "exit" ||
                   fun->getName().str() == "__assert_rtn" ||
                   fun->getName().str() == "__assert_fail" );
}

inline bool isProgExitCall(const CallSite cs) {
    return isProgExitFunction(getCallee(cs));
}

inline bool isProgExitCall(const Instruction *inst) {
    return isProgExitFunction(getCallee(inst));
}
//@}

/// Function does not have any possible caller in the call graph
//@{
/// Return true if the function does not have a caller (either it is a main function or a dead function)
inline bool isNoCallerFunction (const Function * fun) {
    return isDeadFunction(fun) || isProgEntryFunction(fun);
}

/// Return true if the argument in a function does not have a caller
inline bool ArgInNoCallerFunction (const Value * val) {
    return SVFUtil::isa<Argument>(val)
           && isNoCallerFunction(SVFUtil::cast<Argument>(val)->getParent());
}
//@}

/// Return true if the function has a return instruction reachable from function entry
bool functionDoesNotRet (const Function * fun);

/// Get reachable basic block from function entry
void getFunReachableBBs (const Function * fun, DominatorTree* dt,std::vector<const BasicBlock*>& bbs);

/// Get function exit basic block
/// FIXME: this back() here is only valid when UnifyFunctionExitNodes pass is invoked
inline const BasicBlock* getFunExitBB(const Function* fun) {
    return &fun->back();
}
/// Strip off the constant casts
const Value * stripConstantCasts(const Value *val);

/// Strip off the all casts
Value *stripAllCasts(Value *val) ;

/// Get the type of the heap allocation
const Type *getTypeOfHeapAlloc(const llvm::Instruction *inst) ;

/// Return corresponding constant expression, otherwise return NULL
//@{
inline const ConstantExpr *isGepConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::GetElementPtr)
            return constExpr;
    }
    return NULL;
}

inline const ConstantExpr *isInt2PtrConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::IntToPtr)
            return constExpr;
    }
    return NULL;
}

inline const ConstantExpr *isPtr2IntConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::PtrToInt)
            return constExpr;
    }
    return NULL;
}

inline const ConstantExpr *isCastConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::BitCast)
            return constExpr;
    }
    return NULL;
}

inline const ConstantExpr *isSelectConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::Select)
            return constExpr;
    }
    return NULL;
}

inline const ConstantExpr *isTruncConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::Trunc ||
                constExpr->getOpcode() == Instruction::FPTrunc ||
                constExpr->getOpcode() == Instruction::ZExt ||
                constExpr->getOpcode() == Instruction::SExt ||
                constExpr->getOpcode() == Instruction::FPExt)
            return constExpr;
    }
    return NULL;
}

inline const ConstantExpr *isCmpConstantExpr(const Value *val) {
    if(const ConstantExpr* constExpr = SVFUtil::dyn_cast<ConstantExpr>(val)) {
        if(constExpr->getOpcode() == Instruction::ICmp || constExpr->getOpcode() == Instruction::FCmp)
            return constExpr;
    }
    return NULL;
}
//@}

/// Get the next instructions following control flow
void getNextInsts(const Instruction* curInst, std::vector<const Instruction*>& instList);

/// Get the previous instructions following control flow
void getPrevInsts(const Instruction* curInst, std::vector<const Instruction*>& instList);

/// Get basic block successor position
u32_t getBBSuccessorPos(const BasicBlock *BB, const BasicBlock *Succ);
/// Get num of BB's successors
u32_t getBBSuccessorNum(const BasicBlock *BB);
/// Get basic block predecessor positin
u32_t getBBPredecessorPos(const BasicBlock *BB, const BasicBlock *Pred);
/// Get num of BB's predecessors
u32_t getBBPredecessorNum(const BasicBlock *BB);

/// Return source code including line number and file name from debug information
//@{
std::string  getSourceLoc(const Value *val);
std::string  getSourceLocOfFunction(const Function *F);
//@}

/// Dump sparse bitvector set
void dumpSet(NodeBS To, raw_ostream & O = SVFUtil::outs());

/// Dump points-to set
void dumpPointsToSet(unsigned node, NodeBS To) ;

/// Dump alias set
void dumpAliasSet(unsigned node, NodeBS To) ;

/// Print successful message by converting a string into green string output
std::string sucMsg(std::string msg);

/// Print warning message by converting a string into yellow string output
void wrnMsg(std::string msg);

/// Print error message by converting a string into red string output
//@{
std::string  errMsg(std::string msg);
std::string  bugMsg1(std::string msg);
std::string  bugMsg2(std::string msg);
std::string  bugMsg3(std::string msg);
//@}

/// Print each pass/phase message by converting a string into blue string output
std::string  pasMsg(std::string msg);

/// Print memory usage in KB.
void reportMemoryUsageKB(const std::string& infor, raw_ostream & O = SVFUtil::outs());

/// Get memory usage from system file. Return TRUE if succeed.
bool getMemoryUsageKB(u32_t* vmrss_kb, u32_t* vmsize_kb);

/// Increase the stack size limit
void increaseStackSize();

/// Check whether a file is an LLVM IR file
bool isIRFile(const std::string &filename);

/// Parse argument for multi-module analysis
void processArguments(int argc, char **argv, int &arg_num, char **arg_value,
                      std::vector<std::string> &moduleNameVec);
/*!
 * Compare two PointsTo according to their size and points-to elements.
 * 1. PointsTo with smaller size is smaller than the other;
 * 2. If the sizes are equal, comparing the points-to targets.
 */
inline bool cmpPts (const PointsTo& lpts,const PointsTo& rpts) {
    if (lpts.count() != rpts.count())
        return (lpts.count() < rpts.count());
    else {
        PointsTo::iterator bit = lpts.begin(), eit = lpts.end();
        PointsTo::iterator rbit = rpts.begin(), reit = rpts.end();
        for (; bit != eit && rbit != reit; bit++, rbit++) {
            if (*bit != *rbit)
                return (*bit < *rbit);
        }

        return false;
    }
}

}

//===- SVFModule.h -- SVFModule class-----------------------------------------//
//
//                     SVF: Static Value-Flow Analysis
//
// Copyright (C) <2013-2017>  <Yulei Sui>
//

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//===----------------------------------------------------------------------===//

/*
 * SVFModule.h
 *
 *  Created on: Aug 4, 2017
 *      Author: Xiaokang Fan
 */



class LLVMModuleSet {
public:
    typedef std::vector<Function*> FunctionSetType;
    typedef std::vector<GlobalVariable*> GlobalSetType;
    typedef std::vector<GlobalAlias*> AliasSetType;

    typedef std::map<const Function*, Function*> FunDeclToDefMapTy;
    typedef std::map<const Function*, FunctionSetType> FunDefToDeclsMapTy;
    typedef std::map<const GlobalVariable*, GlobalVariable*> GlobalDefToRepMapTy;

    /// Iterators type def
    typedef FunctionSetType::iterator iterator;
    typedef FunctionSetType::const_iterator const_iterator;
    typedef GlobalSetType::iterator global_iterator;
    typedef GlobalSetType::const_iterator const_global_iterator;
    typedef AliasSetType::iterator alias_iterator;
    typedef AliasSetType::const_iterator const_alias_iterator;

private:
    u32_t moduleNum;
    LLVMContext *cxts;
    std::unique_ptr<Module> *modules;

    FunctionSetType FunctionSet;  ///< The Functions in the module
    GlobalSetType GlobalSet;      ///< The Global Variables in the module
    AliasSetType AliasSet;        ///< The Aliases in the module

    /// Function declaration to function definition map
    FunDeclToDefMapTy FunDeclToDefMap;
    /// Function definition to function declaration map
    FunDefToDeclsMapTy FunDefToDeclsMap;
    /// Global definition to a rep definition map
    GlobalDefToRepMapTy GlobalDefToRepMap;

public:
    /// Constructor
    LLVMModuleSet(const std::vector<std::string> &moduleNameVec);
    LLVMModuleSet(Module *mod);
    LLVMModuleSet(Module &mod);
    LLVMModuleSet() {}

    void build(const std::vector<std::string> &moduleNameVec);

    u32_t getModuleNum() const {
        return moduleNum;
    }

    Module *getModule(u32_t idx) const {
        assert(idx < moduleNum && "Out of range.");
        return modules[idx].get();
    }

    Module &getModuleRef(u32_t idx) const {
        assert(idx < moduleNum && "Out of range.");
        return *(modules[idx].get());
    }

    // Dump modules to files
    void dumpModulesToFile(const std::string suffix);

    /// Fun decl --> def
    bool hasDefinition(const Function *fun) const {
        assert(fun->isDeclaration() && "not a function declaration?");
        FunDeclToDefMapTy::const_iterator it = FunDeclToDefMap.find(fun);
        return it != FunDeclToDefMap.end();
    }

    Function *getDefinition(const Function *fun) const {
        assert(fun->isDeclaration() && "not a function declaration?");
        FunDeclToDefMapTy::const_iterator it = FunDeclToDefMap.find(fun);
        assert(it != FunDeclToDefMap.end() && "has no definition?");
        return it->second;
    }

    /// Fun def --> decl
    bool hasDeclaration(const Function *fun) const {
        assert(!fun->isDeclaration() && "not a function definition?");
        FunDefToDeclsMapTy::const_iterator it = FunDefToDeclsMap.find(fun);
        return it != FunDefToDeclsMap.end();
    }

    const FunctionSetType &getDeclaration(const Function *fun) const {
        assert(!fun->isDeclaration() && "not a function definition?");
        FunDefToDeclsMapTy::const_iterator it = FunDefToDeclsMap.find(fun);
        assert(it != FunDefToDeclsMap.end() && "has no declaration?");
        return it->second;
    }

    /// Global to rep
    bool hasGlobalRep(const GlobalVariable *val) const {
        GlobalDefToRepMapTy::const_iterator it = GlobalDefToRepMap.find(val);
        return it != GlobalDefToRepMap.end();
    }

    GlobalVariable *getGlobalRep(const GlobalVariable *val) const {
        GlobalDefToRepMapTy::const_iterator it = GlobalDefToRepMap.find(val);
        assert(it != GlobalDefToRepMap.end() && "has no rep?");
        return it->second;
    }

    /// Iterators
    ///@{
    iterator begin() {
        return FunctionSet.begin();
    }
    const_iterator begin() const {
        return FunctionSet.begin();
    }
    iterator end() {
        return FunctionSet.end();
    }
    const_iterator end() const {
        return FunctionSet.end();
    }

    global_iterator global_begin() {
        return GlobalSet.begin();
    }
    const_global_iterator global_begin() const {
        return GlobalSet.begin();
    }
    global_iterator global_end() {
        return GlobalSet.end();
    }
    const_global_iterator global_end() const {
        return GlobalSet.end();
    }

    alias_iterator alias_begin() {
        return AliasSet.begin();
    }
    const_alias_iterator alias_begin() const {
        return AliasSet.begin();
    }
    alias_iterator alias_end() {
        return AliasSet.end();
    }
    const_alias_iterator alias_end() const {
        return AliasSet.end();
    }
    ///@}

private:
    void loadModules(const std::vector<std::string> &moduleNameVec);
    void addSVFMain();
    void initialize();
    void buildFunToFunMap();
    void buildGlobalDefToRepMap();
};

class SVFModule {
public:
    typedef LLVMModuleSet::FunctionSetType FunctionSetType;
    typedef LLVMModuleSet::GlobalSetType GlobalSetType;
    typedef LLVMModuleSet::AliasSetType AliasSetType;

    typedef LLVMModuleSet::FunDeclToDefMapTy FunDeclToDefMapTy;
    typedef LLVMModuleSet::FunDefToDeclsMapTy FunDefToDeclsMapTy;
    typedef LLVMModuleSet::GlobalDefToRepMapTy GlobalDefToRepMapTy;

    /// Iterators type def
    typedef FunctionSetType::iterator iterator;
    typedef FunctionSetType::const_iterator const_iterator;
    typedef GlobalSetType::iterator global_iterator;
    typedef GlobalSetType::const_iterator const_global_iterator;
    typedef AliasSetType::iterator alias_iterator;
    typedef AliasSetType::const_iterator const_alias_iterator;

private:
    static LLVMModuleSet *llvmModuleSet;
    static std::string pagReadFromTxt;

public:
    /// Constructors
    SVFModule(const std::vector<std::string> &moduleNameVec) {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet(moduleNameVec);
    }
    SVFModule(Module *mod) {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet(mod);
    }
    SVFModule(Module &mod) {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet(mod);
    }
    SVFModule() {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet;
    }

    static inline LLVMModuleSet *getLLVMModuleSet() {
        if (llvmModuleSet == NULL)
            llvmModuleSet = new LLVMModuleSet;
        return llvmModuleSet;
    }

    static inline void setPagFromTXT(std::string txt) {
        pagReadFromTxt = txt;
    }

    static inline std::string pagFileName() {
        return pagReadFromTxt;
    }

    static inline bool pagReadFromTXT() {
    		if(pagReadFromTxt.empty())
    			return false;
    		else
    			return true;
    }

    static void releaseLLVMModuleSet() {
        if (llvmModuleSet)
            delete llvmModuleSet;
        llvmModuleSet = NULL;
    }

    bool empty() const {
        return getModuleNum() == 0;
    }

    /// Methods from LLVMModuleSet
    u32_t getModuleNum() const {
        return llvmModuleSet->getModuleNum();
    }

    Module *getModule(u32_t idx) const {
        return llvmModuleSet->getModule(idx);
    }

    Module &getModuleRef(u32_t idx) const {
        return llvmModuleSet->getModuleRef(idx);
    }

    // Dump modules to files
    void dumpModulesToFile(const std::string suffix) const {
        llvmModuleSet->dumpModulesToFile(suffix);
    }

    /// Fun decl --> def
    bool hasDefinition(const Function *fun) const {
        return llvmModuleSet->hasDefinition(fun);
    }

    Function *getDefinition(const Function *fun) const {
        return llvmModuleSet->getDefinition(fun);
    }

    /// Fun def --> decl
    bool hasDeclaration(const Function *fun) const {
        return llvmModuleSet->hasDeclaration(fun);
    }

    const FunctionSetType &getDeclaration(const Function *fun) const {
        return llvmModuleSet->getDeclaration(fun);
    }

    /// Global to rep
    bool hasGlobalRep(const GlobalVariable *val) const {
        return llvmModuleSet->hasGlobalRep(val);
    }

    GlobalVariable *getGlobalRep(const GlobalVariable *val) const {
        return llvmModuleSet->getGlobalRep(val);
    }

    /// Iterators
    ///@{
    iterator begin() {
        return llvmModuleSet->begin();
    }
    const_iterator begin() const {
        return llvmModuleSet->begin();
    }
    iterator end() {
        return llvmModuleSet->end();
    }
    const_iterator end() const {
        return llvmModuleSet->end();
    }

    Module *getMainLLVMModule() const {
        return llvmModuleSet->getModule(0);
    }

	const std::string& getModuleIdentifier() const {
		if (pagReadFromTxt.empty()) {
			assert(!empty() && "empty LLVM module!!");
			return getMainLLVMModule()->getModuleIdentifier();
		} else {
			return pagReadFromTxt;
		}
	}

    LLVMContext& getContext() const {
        assert(!empty() && "empty LLVM module!!");
        return getMainLLVMModule()->getContext();
    }

    inline Function* getFunction(StringRef name) const {
        Function* fun = NULL;
        for (u32_t i = 0; i < getModuleNum(); ++i) {
            Module *mod = llvmModuleSet->getModule(i);
            fun = mod->getFunction(name);
            if(fun && !fun->isDeclaration()) {
                return fun;
            }
        }
        return fun;
    }

    global_iterator global_begin() {
        return llvmModuleSet->global_begin();
    }
    const_global_iterator global_begin() const {
        return llvmModuleSet->global_begin();
    }
    global_iterator global_end() {
        return llvmModuleSet->global_end();
    }
    const_global_iterator global_end() const {
        return llvmModuleSet->global_end();
    }

    alias_iterator alias_begin() {
        return llvmModuleSet->alias_begin();
    }
    const_alias_iterator alias_begin() const {
        return llvmModuleSet->alias_begin();
    }
    alias_iterator alias_end() {
        return llvmModuleSet->alias_end();
    }
    const_alias_iterator alias_end() const {
        return llvmModuleSet->alias_end();
    }
    ///@}
};


