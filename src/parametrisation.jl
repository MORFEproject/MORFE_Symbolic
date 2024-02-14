"""
Data structure containing all the information about the parametrisation.
Naming conventions are the same as in the paper:
A. Vizzaccaro, G. Gobat, A. Frangi, C. Touze'. Direct parametrisation of invariant manifolds 
for non-autonomous forced systems including superharmonic resonances (2023).
"""
mutable struct parametrisation_struct
    RHS_d::Matrix{Sym}
    RHS_Q::Matrix{Sym}
    W::Matrix{Sym}
    Ws::Matrix{Sym}
    Wr::Matrix{Sym}
    Wmodal::Matrix{Sym}
    f::Matrix{Sym}
    fs::Matrix{Sym}
    fr::Matrix{Sym}
    σ::Matrix{Sym}
    res::Matrix{Sym}
    AYR::Matrix{Sym}
    YLᵀA::Matrix{Sym}
    subs::Vector{Dict{Sym, Sym}}
    n_full::Int
    n_rom::Int
    n_sets::Int
    n_aut::Int
    n_osc::Int
    order::Int
end

"""
Initializes the parametrisation_struct
"""
function init_parametrisation_struct(n_full::Int64,n_rom::Int64,n_sets::Int64,n_aut::Int64,n_osc::Int64,order::Int64)
    mat_Nxns = sympy.zeros(n_full,n_sets);                 # RHS and W    
    mat_nrxns = sympy.zeros(n_rom,n_sets);                 # f
    vec_1xns = sympy.zeros(1,n_sets);                      # σ
    mat_naxns = sympy.zeros(n_aut,n_sets);                 # resonances of index_set with λₐᵤₜ
    mat_Nxna = sympy.zeros(n_full,n_aut);                  # borders of homological matrix (YLᵀ*A and A*YR)
    vec_subs = [Dict(Sym(0)=>Sym(0));Dict(Sym(0)=>Sym(0))] # subs 
    return parametrisation_struct(
        0*mat_Nxns,0*mat_Nxns,0*mat_Nxns,0*mat_Nxns,0*mat_Nxns,0*mat_Nxns,
        0*mat_nrxns,0*mat_nrxns,0*mat_nrxns,
        0*vec_1xns,0*mat_naxns,
        0*mat_Nxna,transpose(mat_Nxna),
        vec_subs,
        n_full,n_rom,n_sets,n_aut,n_osc,order)
end