module MORFE_Symbolic

    using SymPy
    using Combinatorics
    using LinearAlgebra
    using MathLink

    include("right_hand_side.jl")
    include("basic_functionalities.jl")
    include("realification.jl")
    include("matcont.jl")
    
    export create_gen_vec,create_pos_vec
    export mysimp,mysub
    export multiexponent_struct,init_multiexponent_struct
    export system_struct,extract_Lin,extract_Quad,Q_fun
    export parametrisation_struct,init_parametrisation_struct
    export fill_RHS_quad!,fill_RHS_dyn!,fill_RHS_lin!
    export compute_order_zero,generalised_eigenproblem,solve_homological!
    export substitutions!, reduced_dynamics_substitutions!, nonlinear_mappings_substitutions!, Mathematica_output
    export polar_realification, cartesian_realification!, backbone_CNF, physical_amplitudes_CNF
    export modal_coordinates_from_physical_coordinates!
    export matcont 

    #~~~~~~~~~~~~~~~~~#
    #           compute                #
    #~~~~~~~~~~~~~~~~~#

    function compute_order_zero(sys::system_struct,sol::Int64)
        N = size(sys.A)[1]
        if sys.C0 == sympy.zeros(N,1)[:,1]
            x0 = sys.C0
        else
            x = create_gen_vec("x",N)
            x0_sol = solve(Q_fun(x,x,sys.Q)+sys.B*x+sys.C0)
            println("Static Solutions:")
            println(x0_sol)
            x0 = sympy.zeros(N,1)[:,1];
            for in = 1:N
                x0[in] = get(x0_sol[sol],x[in],"-")
            end
        end
        return x0
    end

    function generalised_eigenproblem(sys::system_struct)
        A = sys.A 
        B = sys.B     
        N = size(A)[1] 
        P,D = A.diagonalize() 
        for i = 1:N
            if D[i,i] == 0
                D[i,i] = 1/Sym("∞")
            end
        end
        invA = P*D^-1*P^-1
        invAB = invA*B
        Y,Λ = invAB.diagonalize()   #     A*Y*Λ = B*Y
        println("Eigenproblem calculation right")
        λ = mysimp(diag(Λ))
        println(λ)
        invAᵀ = transpose(P)^-1*D^-1*transpose(P)
        invAᵀBᵀ = invAᵀ*transpose(B)
        X,Λ = invAᵀBᵀ.diagonalize()   #     Aᵀ*X*Λ = Bᵀ*X =>  Λ*Xᵀ*A = Xᵀ*B
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

    function solve_homological!(ind_set,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
        res_ind = DP.res[:,ind_set];
        σs =  Sym("σ"*string(ind_set))#DP.σ[ind_set]
        DP.subs = [DP.subs;Dict(σs=>DP.σ[ind_set])]
        Mat = [   [σs*sys.A-sys.B         DP.AYR*diagm(res_ind)];
                    [   diagm(res_ind)*DP.YLᵀA           diagm(res_ind.-Sym(1))]];
        Vec = [   0*(DP.RHS_d[:,ind_set]+DP.RHS_Q[:,ind_set]);
                    sympy.zeros(DP.n_aut,1)[:,1]];
        for i_full = 1:DP.n_full
            if DP.RHS_d[i_full,ind_set]+DP.RHS_Q[i_full,ind_set] !=0
                Vec[i_full] =  Sym("RHS"*string(i_full)*"0"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("RHS"*string(i_full)*"0"*string(ind_set))=>DP.RHS_d[i_full,ind_set]+DP.RHS_Q[i_full,ind_set])]
            end
        end
        Sol = Mat\Vec
        DP.W[:,ind_set] = Sol[1:DP.n_full]
        for i_full=1:DP.n_full
            if DP.W[i_full,ind_set] != 0
                DP.Ws[i_full,ind_set] = Sym("W"*string(i_full)*"0"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"0"*string(ind_set))=>DP.W[i_full,ind_set])]
            end
        end
        DP.f[1:DP.n_aut,ind_set] = Sol[DP.n_full+1:DP.n_full+DP.n_aut]
        for i_aut=1:DP.n_aut
            if DP.f[i_aut,ind_set] != 0
                DP.fs[i_aut,ind_set] = Sym("f"*string(i_aut)*"0"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("f"*string(i_aut)*"0"*string(ind_set))=>DP.f[i_aut,ind_set])]
            end
        end
    end

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

    function Mathematica_output(DP::parametrisation_struct, aexp::multiexponent_struct, directory, file_basename;
        print_reduced_dynamics = false, print_nonlinear_mappings = false, 
        print_cartesian_realified_reduced_dynamics = false, print_cartesian_realified_nonlinear_mappings = false,
        print_polar_realified_reduced_dynamics = false, real = nothing, imaginary = nothing,
        print_backbone = false, omega_rho = nothing, print_damping = false, xi_rho = nothing,
        print_physical_amplitudes = false, ampli_rho = nothing)
        mathematica_code = sympy.printing.mathematica.mathematica_code

        println("Mathematica output started")

        isdir(directory) || mkdir(directory)
        path = joinpath(directory, file_basename * "_variables.nb")
        open(path, "w") do file
            write(file, "naut = $(DP.n_aut);\n")
            write(file, "nfull = $(DP.n_full);\n")
            write(file, "nsets = $(DP.n_sets);\n")
            write(file, "order = $(DP.order);\n")
            for i in eachindex(DP.subs)
                for key in keys(DP.subs[i])
                    if key != 0
                        math_code = mathematica_code(key) * " = FullSimplify[" * mathematica_code(DP.subs[i][key]) * "];\n"
                        write(file, math_code) 
                    end
                end
            end
            write(file, "monoms = ConstantArray[0,$(DP.n_sets)];\n")
            for i_set = 1:DP.n_sets
                z = [Sym("z$(i)") for i=1:DP.n_rom]
                monom = prod(z .^ aexp.mat[:,i_set])
                monom = mathematica_code(monom)
                write(file, "monoms[[$(i_set)]] = " * monom * ";\n")
            end
            if print_reduced_dynamics
                write(file, "f = ConstantArray[0,{$(DP.n_rom),$(DP.n_sets)}];\n")
                for i_var = 1:DP.n_rom
                    for i_set = 1:DP.n_sets
                        write(file, "f[[$(i_var),$(i_set)]] = FullSimplify[" * mathematica_code(DP.f[i_var,i_set]) * "];\n")
                    end
                end
            end
            if print_nonlinear_mappings
                write(file, "W = ConstantArray[0,{$(DP.n_full),$(DP.n_sets)}];\n")
                for i_var = 1:DP.n_full
                    for i_set = 1:DP.n_sets
                        write(file, "W[[$(i_var),$(i_set)]] = FullSimplify[" * mathematica_code(DP.W[i_var,i_set]) * "];\n")
                    end
                end
            end
            if print_cartesian_realified_reduced_dynamics
                write(file, "fr = ConstantArray[0,{$(DP.n_rom),$(DP.n_sets)}];\n")
                for i_var = 1:DP.n_rom
                    for i_set = 1:DP.n_sets
                        write(file, "fr[[$(i_var),$(i_set)]] = FullSimplify[" * mathematica_code(DP.fr[i_var,i_set]) * "];\n")
                    end
                end
            end
            if print_cartesian_realified_nonlinear_mappings
                write(file, "Wr = ConstantArray[0,{$(DP.n_full),$(DP.n_sets)}];\n")
                for i_var = 1:DP.n_full
                    for i_set = 1:DP.n_sets
                        write(file, "Wr[[$(i_var),$(i_set)]] = FullSimplify[" * mathematica_code(DP.Wr[i_var,i_set]) * "];\n")
                    end
                end
            end
            if print_backbone
                write(file, "omegaRho = ConstantArray[0,$(DP.order)];\n")
                for i = 1:DP.order
                    write(file, "omegaRho[[$(i)]] = FullSimplify[" * mathematica_code(omega_rho[i]) * "];\n")
                end
            end
            if print_damping
                write(file, "xiRho = ConstantArray[0,$(DP.order)];\n")
                for i = 1:DP.order
                    write(file, "xiRho[[$(i)]] = FullSimplify[" * mathematica_code(xi_rho[i]) * "];\n")
                end
            end
            if print_physical_amplitudes
                write(file, "ampliRho = ConstantArray[0,$(DP.order)];\n")
                for i = 1:DP.order
                    write(file, "ampliRho[[$(i)]] = FullSimplify[" * mathematica_code(ampli_rho[i]) * "];\n")
                end
            end
        end

        if print_reduced_dynamics
            path = joinpath(directory, file_basename * "_reduced_dynamics.nb")
            open(path, "w") do file
                write(file, "fres = ConstantArray[0, naut];\n")
                write(file, "For[ivar = 1, ivar <= naut, ivar++,\n")
                write(file, "   For[iset = 1, iset <= nsets, iset++,\n")
                write(file, "       fres[[ivar]] += f[[ivar, iset]]*monoms[[iset]];\n")
                write(file, "   ];\n")
                write(file, "];\n")
            end
        end

        if print_nonlinear_mappings
            path = joinpath(directory, file_basename * "_nonlinear_mappings.nb")
            open(path, "w") do file
                write(file, "Wres = ConstantArray[0, nfull];\n")
                write(file, "For[ivar = 1, ivar <= nfull, ivar++,\n")
                write(file, "   For[iset = 1, iset <= nsets, iset++,\n")
                write(file, "       Wres[[ivar]] += W[[ivar, iset]]*monoms[[iset]];\n")
                write(file, "   ];\n")
                write(file, "];\n")
            end
        end

        if print_cartesian_realified_reduced_dynamics
            path = joinpath(directory, file_basename * "_cartesian_realified_reduced_dynamics.nb")
            open(path, "w") do file
                write(file, "frres = ConstantArray[0, naut];\n")
                write(file, "For[ivar = 1, ivar <= naut, ivar++,\n")
                write(file, "   For[iset = 1, iset <= nsets, iset++,\n")
                write(file, "       frres[[ivar]] += fr[[ivar, iset]]*monoms[[iset]];\n")
                write(file, "   ];\n")
                write(file, "];\n")
            end
        end

        if print_cartesian_realified_nonlinear_mappings
            path = joinpath(directory, file_basename * "_cartesian_realified_nonlinear_mappings.nb")
            open(path, "w") do file
                write(file, "Wrres = ConstantArray[0, nfull];\n")
                write(file, "For[ivar = 1, ivar <= nfull, ivar++,\n")
                write(file, "   For[iset = 1, iset <= nsets, iset++,\n")
                write(file, "       Wrres[[ivar]] += Wr[[ivar, iset]]*monoms[[iset]];\n")
                write(file, "   ];\n")
                write(file, "];\n")
            end
        end

        if print_backbone
            path = joinpath(directory, file_basename * "_backbone.nb")
            open(path, "w") do file
                write(file, "backbone = 0;\n")
                write(file, "For[i = 1, i <= order, i++,\n")
                write(file, "   backbone += omegaRho[[i]];\n")
                write(file, "];\n")
            end
        end

        if print_damping
            path = joinpath(directory, file_basename * "_damping.nb")
            open(path, "w") do file
                write(file, "damping = 0;\n")
                write(file, "For[i = 1, i <= order, i++,\n")
                write(file, "   damping += xiRho[[i]];\n")
                write(file, "];\n")
            end
        end

        if print_physical_amplitudes
            path = joinpath(directory, file_basename * "_physical_amplitudes.nb")
            open(path, "w") do file
                write(file, "physicalAmplitudes = 0;\n")
                write(file, "For[i = 1, i <= order, i++,\n")
                write(file, "   physicalAmplitudes += ampliRho[[i]];\n")
                write(file, "];\n")
            end
        end
        
        println("Mathematica output finished")
    end

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




