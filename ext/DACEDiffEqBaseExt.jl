module DACEDiffEqBaseExt
using DACE, DiffEqBase

DiffEqBase.value(a::DACE.DA) = DACE.cons(a)

end
