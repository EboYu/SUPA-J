%module DDAPass_SWIG
%{
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAPass.h"
%}

%import "llvm/Pass.h"
%import "../Util/BasicTypes.i"
%include <cpointer.i>
%include <stdint.i>
%include <std_set.i>
%include <std_vector.i>
%pointer_class(PointerAnalysis,PointerAnalysisP)
%pointer_class(SVFGSCC,SVFGSCCP)
%pointer_class(DDAClient,DDAClientP)
%pointer_class(SVFGEdge,SVFGEdgeP)
%pointer_class(Value,ValueP)
%pointer_class(SVFG,SVFGP)
namespace std {
%template(IntVector) vector<int>;
%template(DoubleVector) vector<double>;
}
/*
 * @file: DDAPass.h
 * @author: Yulei Sui
 * @date: 01/07/2014
 * @version: 1.0
 *
 */


#ifndef WPA_H_
#define WPA_H_

#include "MemoryModel/PointerAnalysis.h"
#include "DDA/DDAClient.h"
#include "Util/SCC.h"
#include "Util/BasicTypes.h"

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


#endif /* WPA_H_ */
