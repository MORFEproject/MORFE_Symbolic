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
    BYR::Matrix{Sym}
    YLᵀB::Matrix{Sym}
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

"""
Function to define the parametrisation and multiexponent data structures. Takes as input arguments:\\
    - λ₀: Array with symbolic variables representing the eigenvalues of the system.\\
    - n_rom: Number of variables in the reduced order model;\\
    - n_full: Number of variables in the full system;\\
    - n_aut: Number of variables in the autonomous part of the equations;\\
    - n_osc: Number of variables in the second order form of the system;\\
    - o: Order of the parametrisation;\\
    - style: String defining the style of parametrisation;\\
    - resonace_conditions: Array of dictionary pairs defining relationships between frequencies through
    λ₀. For example [λ₀[2] =>-λ₀[1]] indicates that the frequency of the first two eigenvalues are the
    negative of one another.\\
Returns the parametrisation structure DP and the multiexponent structure aexp.
"""
function define_parametrisation(λ₀,n_rom,n_full,n_aut,n_osc,o,style,resonance_conditions)
    aexp = init_multiexponent_struct(n_rom,o)
    DP = init_parametrisation_struct(n_full,n_rom,aexp.n_sets,n_aut,n_osc,o)

    σ₀ = transpose(aexp.mat)*λ₀
    if style == "Graph"
        DP.res = DP.res.+1
    else
        for in = 1:n_aut
            if style == "RNF"            
                DP.res[in,:] = convert(Vector{Int},mysub(abs.(σ₀).-abs(λ₀[in]),resonance_conditions) .== 0)
            end
            if style == "CNF"
                DP.res[in,:] = convert(Vector{Int},mysub(σ₀.-λ₀[in],resonance_conditions) .== 0)
            end
        end
    end

    return DP, aexp
end

"""
Computes the order zero of the parametrisation. Takes as input arguments:\\
    - DP: The parametrisation structure;\\
    - aexp: The multiexponent structure;\\
    - sys: The system of equations structure;\\
Does not return anything, but changes DP.
"""
function compute_order_0_parametrisation!(DP, aexp, sys)
    p0 = 0;                                       # Order of the parametrisation
    n0 = 1;                                       # Number of sets in order 0
    ind_set0 = 1                                  # Local index related to the set of order 0
    ind_setG0 = aexp.get(aexp.get([p0 ind_set0])) # Global index related to the set of order 0

    Y0 = compute_order_zero(sys,1)                # Computes static solution. Read function description for
                                                  # more details.
    DP.W[:,ind_setG0] = Y0
    for i_full=1:DP.n_full
        if DP.W[i_full,ind_setG0] != 0
            DP.Ws[i_full,ind_setG0] = Sym("W"*string(i_full)*"X"*string(ind_setG0))
            DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"X"*string(ind_setG0))=>DP.W[i_full,ind_setG0])]
        end
    end

    # the matrix B is now updated by adding to it 
    # the gradient of the quadratic RHS calculated at Y0
    ∇Q0_fun(Y) = Q_fun(Y,Y0,sys.Q)
    ∇Q0 = extract_Lin(∇Q0_fun,DP.n_full)
    sys.B += ∇Q0
end

"""
Computes the order one of the parametrisation. Takes as input arguments:\\
    - DP: The parametrisation structure;\\
    - aexp: The multiexponent structure;\\
    - sys: The system of equations structure;\\
    - eigenvalues_order: Variable defining the order of the eigenvalues. Needs to be of a type
    such that indexation of an array with it is possible (can be an integer, Vector, StepRange, etc);\\
    - nonautonomous_eigenvalues: Column vector with the expressions for the eigenvalues of the
    nonautonomous part.\\
Does not return anything, but changes DP.
"""
function compute_order_1_parametrisation!(DP, aexp, sys, eigenvalues_order, nonautonomous_eigenvalues = nothing)
    Λ,YR,YL = generalised_eigenproblem(sys)
    
    # Master coordinates
    λ = Λ[eigenvalues_order]
    yR = YR[:,eigenvalues_order]
    yL = YL[:,eigenvalues_order]
    # Any choice is possible but the sorting is not known before launching the script!
    # One must then look at the eigenvalues sorting, then choose the masters
    
    # Right eigenvalues normalisation
    for i = 1:size(yR)[2]
        for j = 1:size(yR)[1]
            if abs(yR[j,i]) > 0
                yR[:,i] = yR[:,i]/(yR[j,i])
                break
            end
        end
    end

    # Filling the parametrisation
    p1 = 1;

    # Autonomous part
    for ind_set1 = 1:DP.n_aut
        ind_setG1 = aexp.get(aexp.get([p1 ind_set1]))
        DP.W[:,ind_setG1] = yR[:,ind_set1]
        for i_full=1:DP.n_full
            if DP.W[i_full,ind_setG1] != 0
                DP.Ws[i_full,ind_setG1] = Sym("W"*string(i_full)*"X"*string(ind_setG1))
                DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"X"*string(ind_setG1))=>DP.W[i_full,ind_setG1])]
            end
        end
        DP.f[ind_set1,ind_setG1] = λ[ind_set1]
        DP.fs[ind_set1,ind_setG1] = Sym("λ"*string(ind_set1))
        DP.subs = [DP.subs;Dict(Sym("λ"*string(ind_set1))=>λ[ind_set1])]
        # compute the matrix B*yR[aut]
        # which will be used for the top right border of the homological matrix
        yRs = 0*yR[:,ind_set1]
        for i_full=1:DP.n_full
            if yR[i_full,ind_set1] != 0
                yRs[i_full] = Sym("yR"*string(i_full)*"X"*string(ind_set1))
                DP.subs = [DP.subs;Dict(Sym("yR"*string(i_full)*"X"*string(ind_set1))=>yR[i_full,ind_set1])]
            end
        end
        DP.BYR[:,ind_set1] = sys.B*yRs
        # compute the matrix yL[aut]ᵀ*B
        # which will be used for the bottom left border of the homological matrix
        yLsᵀ = 0*transpose(yL[:,ind_set1])
        for i_full=1:DP.n_full
            if yL[i_full,ind_set1] != 0
                yLsᵀ[i_full] = Sym("yL"*string(i_full)*"X"*string(ind_set1))
                DP.subs = [DP.subs;Dict(Sym("yL"*string(i_full)*"X"*string(ind_set1))=>yL[i_full,ind_set1])]
            end
        end
        DP.YLᵀB[ind_set1,:] = yLsᵀ*sys.B 
    end

    # Nonautonomous part
    n_nonaut = DP.n_rom - DP.n_aut
    if n_nonaut == 0
        DP.σ = transpose(λ)*aexp.mat
    else
        # Augment λ with eigenvalues of the nonautonomous part:
        if nonautonomous_eigenvalues !== nothing
            λ = [λ;nonautonomous_eigenvalues]
        end
        λ = reshape(λ,1,DP.n_rom)
        # Assign the eigenvalues of the nonautonomous part to f:
        DP.f[DP.n_aut+1,aexp.get(aexp.get([p1 DP.n_aut+1]))] = λ[DP.n_aut+1]
        DP.fs[DP.n_aut+1,aexp.get(aexp.get([p1 DP.n_aut+1]))] = Sym("λ"*string(DP.n_aut+1))
        DP.subs = [DP.subs;Dict(Sym("λ"*string(DP.n_aut+1))=>λ[DP.n_aut+1])]
        DP.f[DP.n_aut+2,aexp.get(aexp.get([p1 DP.n_aut+2]))] = λ[DP.n_aut+2]
        DP.fs[DP.n_aut+2,aexp.get(aexp.get([p1 DP.n_aut+2]))] = Sym("λ"*string(DP.n_aut+2))
        DP.subs = [DP.subs;Dict(Sym("λ"*string(DP.n_aut+2))=>λ[DP.n_aut+2])]
        # Assign the C⁺ₑₓₜ and C⁻ₑₓₜ to the RHS of the nonautonomous homological:
        DP.RHS_d[1:DP.n_full,aexp.get(aexp.get([p1 DP.n_aut+1]))] = sys.C⁺ₑₓₜ
        DP.RHS_d[1:DP.n_full,aexp.get(aexp.get([p1 DP.n_aut+2]))] = sys.C⁻ₑₓₜ
        # Calculate σ and assign it to DP:
        DP.σ = λ*aexp.mat
        # σ is equal to:
        # [0;                                             # order 0
        #  λ₁; λ₂; .. λₙᵣₒₘ;                              # order 1
        #  2λ₁; λ₁+λ₂; .. 2λₙᵣₒₘ;                         # order 2
        #  3λ₁; 2λ₁+λ₂; .. 3λₙᵣₒₘ;                        # order 3
        #  ... ]                                          # and so on
        # Solve homological for the linear nonautonomous part:
        for ind_set1_nonaut = DP.n_aut+1:DP.n_aut+n_nonaut
            ind_setG1_nonaut = aexp.get(aexp.get([p1 ind_set1_nonaut]))
            println("solving order "*string(p1)*" and set "*string(ind_setG1_nonaut)*" with exponents:")
            println(aexp.mat[:,ind_setG1_nonaut])
            solve_homological!(ind_setG1_nonaut,DP,sys)
        end
    end
end

"""
Computes the order p of the parametrisation. Takes as input arguments:\\
    - DP: The parametrisation structure;\\
    - aexp: The multiexponent structure;\\
    - sys: The system of equations structure;\\
    - p: The current order.\\
Does not return anything, but changes DP.
"""
function compute_order_p_parametrisation!(DP, aexp, sys, p)
    println("solving order "*string(p))
    fill_RHS_dyn!(p,DP,aexp,sys)
    fill_RHS_quad!(p,DP,aexp,sys)
    np = aexp.get(p)
    for ind_setp = 1:np
        ind_setGp = aexp.get(aexp.get([p ind_setp]))
        println("solving set "*string(ind_setGp)*" with exponents:")
        println(aexp.mat[:,ind_setGp])
        fill_RHS_lin!(aexp.mat[:,ind_setGp],DP,aexp,sys)
        solve_homological!(ind_setGp,DP,sys)
    end
end