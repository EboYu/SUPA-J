%module DDAPass_JAVA
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
%include <inttypes.i>
%include <attribute.i>
%include <exception.i>
%include <cdata.i>
%include <cmalloc.i>
%include <constraints.i>
%include <cwstring.i>
%include <intrusive_ptr.i>
%include <math.i>
%include <std_except.i>
%include <swigarch.i>
%include <swigrun.i>
%include <wchar.i>
%include <shared_ptr.i>

%{
#include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAPass.h"
%}


%import "/home/yinbo/LLVM/llvm-9.0.0.src/include/llvm/Pass.h"
%import "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/BasicTypes.h"
%import "/home/yinbo/disk/workspace/SUPA-J/SVF/include/MemoryModel/PointerAnalysis.h"
%import "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAClient.h"
%import "/home/yinbo/disk/workspace/SUPA-J/SVF/include/Util/SCC.h"
// %pointer_functions(Value, ValueP)
// %pointer_functions(PointerAnalysis,PointerAnalysisP)
// %pointer_functions(SVFG,SVFGP)
// %pointer_functions(SVFGSCC,SVFGSCCP)
// %pointer_functions(SVFGEdge,SVFGEdgeP)
// %pointer_functions(DDAClient,DDAClientP)
%include "/home/yinbo/disk/workspace/SUPA-J/SVF/include/DDA/DDAPass.h"


