%module SUPA

%{
#include "DDA/DDAPass.h"
%}
%include <cpointer.i>
%include <carrays.i>
%pointer_functions(Value, ValueP)
%pointer_functions(PointerAnalysis,PointerAnalysisP)
%pointer_functions(SVFG,SVFGP)
%pointer_functions(SVFGSCC,SVFGSCCP)
%pointer_functions(SVFGEdge,SVFGEdgeP)
%pointer_functions(DDAClient,DDAClientP)
%import "llvm/Pass.h"
%import "../SVF/include/Util/BasicTypes.h"
%include "../SVF/include/DDA/DDAPass.h"
