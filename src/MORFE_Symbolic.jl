module MORFE_Symbolic

    using SymPy
    using Combinatorics
    using LinearAlgebra

    include("basic_functionalities.jl")
    include("system_of_equations.jl")
    include("parametrisation.jl")
    include("multiexponent.jl")
    include("output.jl")
    include("realification.jl")
    include("right_hand_side.jl")
    
    export create_gen_vec,create_real_vec,create_pos_vec
    export mysimp,mysub
    export multiexponent_struct,init_multiexponent_struct
    export system_struct,extract_Lin,extract_Quad,Q_fun
    export parametrisation_struct,init_parametrisation_struct
    export fill_RHS_quad!,fill_RHS_dyn!,fill_RHS_lin!
    export compute_order_zero,generalised_eigenproblem,solve_homological!
    export substitutions!, reduced_dynamics_substitutions!, nonlinear_mappings_substitutions!, Mathematica_output
    export polar_realification, cartesian_realification!, backbone_CNF, physical_amplitudes_CNF
    export modal_coordinates_from_physical_coordinates!
    export reduced_dynamics_latex_output, nonlinear_mappings_latex_output, backbone_output, nonlinear_damping_output
    export physical_amplitudes_output, Mathematica_output
    export partitions_two
    export define_second_order_matrices, define_system
    export define_parametrisation, compute_order_0_parametrisation!
    export compute_order_1_parametrisation!, compute_order_p_parametrisation!

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #             compute              #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    """
    Function to calculate static solution. Finds Y0 such that A.Y0 + (Q.Y0).Y0 + C0 = 0
    ATTENTION: there might be multiple solutions to this problem.
    The second entry of the function selects which solution to extract.
    If C0 = 0, the function compute_order_zero returns automatically Y0 = 0.
    If the function struggles and if the right solution to the equation is known by the user,
    it can be directly input as: DP.W[:,ind_setG0] = [W0₁ ...]
    """
    function compute_order_zero(sys::system_struct,sol::Int64)
        N = size(sys.B)[1]
        if sys.C0 == sympy.zeros(N,1)[:,1]
            x0 = sys.C0
        else
            x = create_gen_vec("x",N)
            x0_sol = solve(Q_fun(x,x,sys.Q)+sys.A*x+sys.C0)
            println("Static Solutions:")
            println(x0_sol)
            x0 = sympy.zeros(N,1)[:,1];
            for in = 1:N
                x0[in] = get(x0_sol[sol],x[in],"-")
            end
        end
        return x0
    end

    """
    Calculates the eigenvalues λ and left and right eigenvalues X and Y
    for the generalized eigenvalue problem defined by (B - λA)Y = X*(B - λA) = 0.
    λ,X,Y fulfil the conditions:
    (1)    sys.B*Y[:,k]*λ[k] = sys.A*Y[:,k]
    (2)    λ[k]*Xᵀ[:,k]*sys.B = Xᵀ[:,k]*sys.A
    For the degrees of freedom of the full system relating to the algebraic part of the DAE,
    the eigenvalues are set to the symbolic variable ∞. The conditions (1) and (2) might not
    be respected for Λ[k] = ∞.
    Normally the eigenvalues are organised as the first n_aux eigenvalues Λ[1:n_aux] are 
    equal to ∞, then the meaningful eigenvalues starts usually sorted by lower ω to higher.
    """
    function generalised_eigenproblem(sys::system_struct)
        A = sys.A 
        B = sys.B     
        N = size(B)[1] 
        P,D = B.diagonalize() 
        for i = 1:N
            if D[i,i] == 0
                D[i,i] = 1/Sym("∞")
            end
        end
        invB = P*D^-1*P^-1
        invBA = invB*A
        Y,Λ = invBA.diagonalize()   #     A*Y*Λ = B*Y
        println("Eigenproblem calculation right")
        λ = mysimp(diag(Λ))
        println(λ)
        invBᵀ = transpose(P)^-1*D^-1*transpose(P)
        invBᵀAᵀ = invBᵀ*transpose(A)
        X,Λ = invBᵀAᵀ.diagonalize()   #     Aᵀ*X*Λ = Bᵀ*X =>  Λ*Xᵀ*A = Xᵀ*B
        println("Eigenproblem calculation left")
        λ = mysimp(diag(Λ))
        println(λ)    
        println("ATTENTION: check that the left and right are in the same order")
        #= for i = 1:N
            a_i = transpose(X[:,i])*A*Y[:,i]  # modal mass
            Y[:,i]/ = sympy.sqrt(a_i)
            X[:,i]/ = sympy.sqrt(a_i)
        end =#
        return [λ,Y,X]
    end

    """
    Solves the homological equation at each order. The results are stored on data structures DP.f and DP.W,
    as a function of symbolic coefficients present at DP.fs and DP.Ws. In order to have results as a function
    of the system parameters, the symbolic coefficients must be substituted using the function substitutions!.
    """
    function solve_homological!(ind_set,DP::parametrisation_struct,sys::system_struct)
        res_ind = DP.res[:,ind_set];
        σs =  Sym("σ"*string(ind_set))#DP.σ[ind_set]
        DP.subs = [DP.subs;Dict(σs=>DP.σ[ind_set])]
        Mat = [   [σs*sys.B-sys.A         DP.BYR*diagm(res_ind)];
                    [   diagm(res_ind)*DP.YLᵀB           diagm(res_ind.-Sym(1))]];
        Vec = [   0*(DP.RHS_d[:,ind_set]+DP.RHS_Q[:,ind_set]);
                    sympy.zeros(DP.n_aut,1)[:,1]];
        for i_full = 1:DP.n_full
            if DP.RHS_d[i_full,ind_set]+DP.RHS_Q[i_full,ind_set] !=0
                Vec[i_full] =  Sym("RHS"*string(i_full)*"X"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("RHS"*string(i_full)*"X"*string(ind_set))=>DP.RHS_d[i_full,ind_set]+DP.RHS_Q[i_full,ind_set])]
            end
        end
        Sol = Mat\Vec
        DP.W[:,ind_set] = Sol[1:DP.n_full]
        for i_full=1:DP.n_full
            if DP.W[i_full,ind_set] != 0
                DP.Ws[i_full,ind_set] = Sym("W"*string(i_full)*"X"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"X"*string(ind_set))=>DP.W[i_full,ind_set])]
            end
        end
        DP.f[1:DP.n_aut,ind_set] = Sol[DP.n_full+1:DP.n_full+DP.n_aut]
        for i_aut=1:DP.n_aut
            if DP.f[i_aut,ind_set] != 0
                DP.fs[i_aut,ind_set] = Sym("f"*string(i_aut)*"X"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("f"*string(i_aut)*"X"*string(ind_set))=>DP.f[i_aut,ind_set])]
            end
        end
    end

    """
    Performs substitutions on vector DP.subs. Should be followed by a call of reduced_dynamics_substitutions!
    and/or nonlinear_mappings_substitutions!.
    """
    function substitutions!(DP::parametrisation_struct,substitutions)
        println("Substituting values:")
        t1 = time_ns()
        for i in eachindex(DP.subs)
            for key in keys(DP.subs[i])
                substituted = mysub([DP.subs[i][key]],DP.subs[1:i-1])
                for substitute in substitutions
                    substituted = mysub(substituted, substitute)
                end
                DP.subs[i][key] = substituted[1]
            end
            println("   Substitution $(i)/$(length(DP.subs)) OK!")
        end
        t2 = time_ns()
        println("Elapsed time: $((t2-t1)/1.0e9) s")
        println("")
    end

    """
    Performs substitutions on the reduced dynamics coefficients. Should be called after substitutions!.
    """
    function reduced_dynamics_substitutions!(DP::parametrisation_struct,substitutions)
        println("Substituting reduced dynamics:")
        t1 = time_ns()
        for i_ord=1:length(DP.f[1,:])
            for i_var=1:DP.n_rom
                substituted = mysub([DP.f[i_var,i_ord]],DP.subs[end:-1:1])
                for substitute in substitutions
                    substituted = mysub(substituted, substitute)
                end
                DP.f[i_var,i_ord] = substituted[1]
            end
            println("   Set $(i_ord)/$(length(DP.f[1,:])) OK!")
        end
        t2 = time_ns()
        println("Elapsed time: $((t2-t1)/1.0e9) s")
        println("")
    end
    
    """
    Performs substitutions on the nonlinear mapping coefficients. Should be called after substitutions!.
    """
    function nonlinear_mappings_substitutions!(DP::parametrisation_struct,substitutions)
        println("Substituting nonlinear mappings:")
        t1 = time_ns()
        for i_ord=1:length(DP.W[1,:])
            for i_var=1:DP.n_full
                substituted = mysub([DP.W[i_var,i_ord]],DP.subs[end:-1:1])
                for substitute in substitutions
                    substituted = mysub(substituted, substitute)
                end
                DP.W[i_var,i_ord] = substituted[1]
            end
            println("   Set $(i_ord)/$(length(DP.f[1,:])) OK!")
        end
        t2 = time_ns()
        println("Elapsed time: $((t2-t1)/1.0e9) s")
        println("")
    end
end




