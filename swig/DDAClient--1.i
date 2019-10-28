%module DDAClient

%{
#include "DDA/DDAClient.h"
%}

%include <cpointer.i>
%pointer_functions(Value, ValueP)
%pointer_functions(PointerAnalysis,PointerAnalysisP)
%pointer_functions(SVFG,SVFGP)
%pointer_functions(SVFGSCC,SVFGSCCP)
%pointer_functions(SVFGEdge,SVFGEdgeP)
%pointer_functions(DDAClient,DDAClientP)

%include "llvm/Pass.h"
%include "../SVF/include/Util/BasicTypes.h"
%include "../SVF/include/DDA/DDAClient.h"
