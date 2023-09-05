module MORFE_Symbolic

    using SymPy
    using Combinatorics
    using LinearAlgebra
    using Latexify
    
    export create_gen_vec,create_pos_vec
    export mysimp,mysub
    export multiexponent_struct,init_multiexponent_struct
    export system_struct,extract_Lin,extract_Quad,Q_fun
    export parametrisation_struct,init_parametrisation_struct
    export fill_RHS_quad!,fill_RHS_dyn!,fill_RHS_lin!
    export compute_order_zero,generalised_eigenproblem,solve_homological!    
    export reduced_dynamics_latex_output, nonlinear_mappings_latex_output
    
    #~~~~~~~~~~~~~~~~~#
    #           generic vectors       #
    #~~~~~~~~~~~~~~~~~#

    function create_gen_vec(str,n)
        z = sympy.zeros(n,1)[:,1]
        for i = 1:n
            z[i] = symbols(str*'_'*string(i))
        end
        return z
    end

    function create_pos_vec(str,n)
        z = sympy.zeros(n,1)[:,1]
        for i = 1:n
            z[i] = symbols(str*'_'*string(i),positive = true)
        end
        return z
    end

    #~~~~~~~~~~~~~~~~~#
    #           manipulation          #
    #~~~~~~~~~~~~~~~~~#

    function mysimp(input)
        exp = reshape(input,prod(size(input)))
        for i = 1:length(exp)
            exp[i] = simplify(exp[i])
        end
        return reshape(exp,size(input))
    end

    function mysub(input,sub)
        exp = reshape(input,prod(size(input)))
        for i = 1:length(exp)
            for j = 1:length(sub)
                exp[i] = exp[i] |> subs(sub[j])
            end
            for j = 1:length(sub)
                exp[i] = expand(exp[i]) |> subs(sub[j])
                exp[i] = simplify(exp[i])
            end
        end
        return reshape(exp,size(input))
    end

    #~~~~~~~~~~~~~~~~~#
    #           multiexponent        #
    #~~~~~~~~~~~~~~~~~#

    function create_multiexp(n::Int,o::Int)
        #= 
        a is a dictionary containing both forward and inverse mappings 
        from a certain multiexponent set [a1;a2;a3..;ap] 
        to its "position" i and order p
        but also its opposite: from [p i] to the set
        moreover, it also contains the number of sets for a certain p

        get_a is a function that interrogates the dictionary
        if given a vector [a1;a2;..;ap] will give the 1x2 matrix [p i]
        if given a 1x2 matrix [p i] will give the vector [a1;a2;..;ap]
        if given an integer p will give the number of vectors [a1;a2;..;ap] of that order 
        =#
        mltexp = collect(multiexponents(n,0));
        MltExp = Matrix{Int}(undef, n,1);
        for in = 1:n;        MltExp[in,1] = mltexp[1][in];    end
        a = Dict(mltexp[1] => 1)
        a = merge(a, Dict([0 1] =>mltexp[1]))
        a = merge(a,Dict(0 =>1))
        n_sets = 1
        for p = 1:o
            mltexp = collect(multiexponents(n,p));
            mltexp_mat = Matrix{Int}(undef, n,length(mltexp))
            for in = 1:n;  for im = 1:length(mltexp);    mltexp_mat[in,im] = mltexp[im][in];    end;end
            MltExp = [MltExp mltexp_mat];
            a = merge(a,Dict(MltExp[:,n_sets+i] => n_sets+i for i = 1:length(mltexp)))
            a = merge(a,Dict([p i] => MltExp[:,n_sets+i] for i = 1:length(mltexp)))
            a = merge(a,Dict(p =>length(mltexp)))
            n_sets = n_sets+length(mltexp);
        end
        function get_a(a_set)
            return get(a,a_set,"-")
        end
        return get_a,a,MltExp,n_sets
    end

    function partitions_two(p::Int)
        # partitions gives a list of all the 2 integers whose sum is p
        # except the trivial one [p,0] and [0,p]
        #
        #
        # the part vector will be for example:
        #= 
        
        partitions_two(5)    = 
        4-element Vector{Vector{Int64}}:
        [4, 1]
        [3, 2]
        [2, 3]
        [1, 4]
        
        =#
        part = collect(partitions(p,2))
        if iseven(p)
            for k = length(part)-1:-1:1
                part = [part;[part[k][end:-1:1]]]
            end
        else
            for k = length(part):-1:1
                part = [part;[part[k][end:-1:1]]]
            end
        end
        return part
    end

    """
    Struct to manage multiexponents. It is constituted of the following:
        - dic: A dictionary containing pairs of (key, value). There are three kinds of keys:
            - [deg local_index]      : Used to access the monomial of degree deg with number local_index.
                                       For example [4, 3] accesses the 3rd monomial of degree 4.
                                       The returned value is an array, in the format of next item in this list.
            - [deg_1, ..., deg_n_rom]: Used to access of the monomial z_1^deg_1 * ... * z_n^deg_n_rom,
                                       where n_rom represents the number of master coordinates.
                                       The returned value is an integer, the global index of the monomial.
            - deg                    : Returns an integer, the number of multiexponents of degree deg.
        - get: A function to get values from dic.
        - mat: An n_rom x n_sets matrix, in which each column is a multiexponent.
        - n_sets: The number of monoms with n_rom variables up to order o.
    """
    mutable struct multiexponent_struct
        get
        dic::Dict{Any,Any}
        mat::Matrix{Int}
        n_sets::Int
    end

    """
    Initializes the multiexponent struct for a model with n_rom master coordinates and
    parametrization up to order o.
    """
    function init_multiexponent_struct(n_rom::Int,o::Int)
        g,d,v,n_sets = create_multiexp(n_rom,o)
        return multiexponent_struct(g,d,v,n_sets)
    end

    #~~~~~~~~~~~~~~~~~#
    #           system original       #
    #~~~~~~~~~~~~~~~~~#

    mutable struct system_struct
        A::Matrix{Sym}
        B::Matrix{Sym}
        Q::Vector{Vector{Sym}}
        C0::Vector{Sym}
        C⁺ₑₓₜ::Vector{Sym}
        C⁻ₑₓₜ::Vector{Sym}
    end

    function extract_Quad(QuadFun,N::Int64)
        Q_tens = [[Sym(0),0,0,0]]
        Y = create_gen_vec("Y",N)
        F = QuadFun(Y)
        for i = 1:N
            for j = 1:N
                for k = 1:N
                    gⁱⱼₖ = diff(diff( F[i] ,Y[j]),Y[k])/Sym(2);
                    if gⁱⱼₖ != 0
                        push!(Q_tens,[i,j,k,gⁱⱼₖ])
                    end
                end
            end
        end
        return Q_tens[2:end]
    end

    function extract_Lin(LinFun,N::Int64)
        Y = create_gen_vec("Y",N)
        F = LinFun(Y)
        L_mat = sympy.zeros(N)
        for i = 1:N
            for j = 1:N
                L_mat[i,j] = diff( F[i] ,Y[j]);
            end
        end
        return L_mat
    end

    function Q_fun(Y1::Vector{Sym},Y2::Vector{Sym},Q::Vector{Vector{Sym}})
        F = Y1*0
        for entry = 1:length(Q)
            i = convert(Int,Q[entry][1])
            j = convert(Int,Q[entry][2])
            k = convert(Int,Q[entry][3])
            F[i]+= Q[entry][4]*Y1[j]*Y2[k]
        end
        return F
    end

    #~~~~~~~~~~~~~~~~~#
    #           parametrisation      #
    #~~~~~~~~~~~~~~~~~#

    mutable struct parametrisation_struct
        RHS_d::Matrix{Sym}
        RHS_Q::Matrix{Sym}
        W::Matrix{Sym}    
        Ws::Matrix{Sym}
        fs::Matrix{Sym}
        f::Matrix{Sym}
        σ::Matrix{Sym}
        res::Matrix{Sym}
        AYR::Matrix{Sym}
        YLᵀA::Matrix{Sym}
        subs::Vector{Dict{Sym, Sym}}
        n_full::Int
        n_rom::Int
        n_sets::Int
        n_aut::Int
    end

    function init_parametrisation_struct(n_full::Int64,n_rom::Int64,n_sets::Int64,n_aut::Int)
        mat_Nxns = sympy.zeros(n_full,n_sets); # RHS and W    
        mat_nrxns = sympy.zeros(n_rom,n_sets);# f
        vec_1xns = sympy.zeros(1,n_sets);          # σ
        mat_naxns = sympy.zeros(n_aut,n_sets);   # resonances of index_set with λₐᵤₜ
        mat_Nxna = sympy.zeros(n_full,n_aut);    # borders of homological matrix (YLᵀ*A and A*YR)
        vec_subs=[Dict(Sym(0)=>Sym(0));Dict(Sym(0)=>Sym(0))]
        #
        return parametrisation_struct(
            0*mat_Nxns,0*mat_Nxns,0*mat_Nxns,0*mat_Nxns,
            0*mat_nrxns,0*mat_nrxns,
            0*vec_1xns,0*mat_naxns,
            0*mat_Nxna,transpose(mat_Nxna),
            vec_subs,
            n_full,n_rom,n_sets,n_aut)
    end

    #~~~~~~~~~~~~~~~~~#
    #           RHS                       #
    #~~~~~~~~~~~~~~~~~#

    function fill_RHS_quad!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
        part = partitions_two(p)
        # partitions gives a list of all the 2 integers whose sum is p
        # expect the trivial one [p,0] and [0,p]
        # the trivial ones must be taken into account in the LHS
        # because they change B to B+∇Q_fun(W[0])
        for k_part = 1:length(part)
            p_L = part[k_part][1];
            p_R = part[k_part][2];
            n_L = aexp.get(p_L);
            n_R = aexp.get(p_R);
            for i_L = 1:n_L
                for i_R = 1:n_R
                    set_L = aexp.get([p_L i_L])
                    set_R = aexp.get([p_R i_R])
                    # the case of symmetric Q is not treated for generality
                    #  RHS_quad(set_L+set_R) += Q(W_map(set_L), W_map(set_R))
                    ind_RHS = aexp.get(set_L+set_R)
                    psiL_ind = aexp.get(set_L)
                    psiR_ind = aexp.get(set_R)
                    DP.RHS_Q[:,ind_RHS] += Q_fun(DP.Ws[:,psiL_ind],DP.Ws[:,psiR_ind],sys.Q)
                end
            end
        end
        return nothing
    end

    function fill_RHS_dyn!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
        part = partitions_two(p+1)    
        # partitions gives a list of all the 2 integers whose sum is p+1
        # expect the trivial one [p+1,0] and [0,p+1]
        # those do not appear because f[0] = 0 and ∂ₜW[0] = 0
        # the outer ones [p,1] and [1,p]
        # are the ones that give:
        # 
        # LHS_dyn(set_current - e_s + e_r) = A*W(set_current)*f(s,r))*set_current[s]
        # in the autonomous case:
        # LHS_dyn(set_current) = A*W(set_current)*f(s,s))*set_current[s]
        # with σ = sum  f(s,s))*set_current[s]
        # 
        # and:
        # LHS_dyn(e_s - e_s + set_current) = A*W(e_s)*f(s,set_current))
        #
        # here I only take the inner ones:
        part = part[2:end-1][:]   
        # NB for forced this has to change! 
        # f^{s}_{s+nA} with  s \in nA is nonzero but ...
        # the p-mapping W is still known from previous calculations of the automonous
        for k_part = 1:length(part)        
            p_L = part[k_part][1];# order of the mapping
            p_R = part[k_part][2];# order of the dynamics
            n_L = aexp.get(p_L);# number of sets in the mappings
            n_R = aexp.get(p_R);# number of sets in the dynamics
            for i_L = 1:n_L
                for i_R = 1:n_R
                    set_L = aexp.get([p_L i_L])# get the set for the mappings
                    set_R = aexp.get([p_R i_R])# get the set for the dynamics
                    for s = 1:length(set_L)       # for s = 1:n_rom
                        if set_L[s]>0               # if the exponent of zₛ is nonzero
                            e_s = 0*set_L;e_s[s] = 1    # create unit vector e_s
                            # the term of the mapping:
                            # W(set_L)*z^(set_L)
                            # when derived in time (along zₛ) gives a term:
                            # W(set_L)*f(s,set_R))*set_L[s]*z^(set_L-e_s+set_R)
                            # then it will go on the RHS with a minus and premult by A:
                            # RHS_dyn(set_L-e_s+set_R) += - A*W(set_L)*f(s,set_R))*set_L[s]
                            ind_f = aexp.get(set_R)
                            fs_set_R = DP.fs[s,ind_f]
                            if fs_set_R != 0    # only fill if the reduced dyn is nonzero
                                ind_RHS = aexp.get(set_L-e_s+set_R)
                                ind_W = aexp.get(set_L)
                                DP.RHS_d[:,ind_RHS] += - sys.A*DP.Ws[:,ind_W]*fs_set_R*set_L[s]
                            end
                        end
                    end                
                end
            end
        end
        return nothing
    end

    function fill_RHS_lin!(set_current::Vector{Int},DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
        # 
        # LHS_dyn(set_dependent - e_s + e_r) = A*W(set_dependent)*f(s,r))*set_dependent[s]
        # sum over s ∈ set_dependent, sum over r ∈ (n_aut+1,n_rom)
        # f(s,r) = δₛᵣ λᵣ    if s and r ∈ (1,n_aut)
        # f(s,r)  != 0      if s ∈ (1,n_aut) and r ∈ (n_aut+1,n_rom)
        # in the autonomous case:
        # LHS_dyn(set_current) = A*W(set_current)*f(s,s))*set_current[s]
        # with σ = sum  f(s,s))*set_current[s]
        # 
        # in the nonautonomous case:
        # for r ∈ n_nonaut
        #   if r ∈ set_current
        #       for s ∈ n_aut
        #           set_dependent = set_current + e_s - e_r   
        #           LHS_dyn(set_current) = A*W(set_dependent)*f(s,r))*set_dependent[s]
        #       end
        #   end
        # end
        for r = 1+DP.n_aut:DP.n_rom
            if set_current[r]>0
                e_r = 0*set_current;e_r[r] = 1    # create unit vector e_r
                ind_f = aexp.get(e_r)
                for s = 1:DP.n_aut                
                    fs_r = DP.f[s,ind_f]                    
                    if fs_r != 0    # only fill if the reduced dyn is nonzero
                        e_s = 0*set_current;e_s[s] = 1    # create unit vector e_s
                        set_dependent = set_current + e_s - e_r
                        ind_W = aexp.get(set_dependent)
                        ind_RHS = aexp.get(set_current)
                        DP.RHS_d[:,ind_RHS] += - sys.A*DP.Ws[:,ind_W]*fs_r*set_dependent[s]
                    end
                end                
            end
        end
        return nothing
    end

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
                Vec[i_full] =  Sym("RHS"*string(i_full)*"|"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("RHS"*string(i_full)*"|"*string(ind_set))=>DP.RHS_d[i_full,ind_set]+DP.RHS_Q[i_full,ind_set])]
            end
        end
        Sol = Mat\Vec
        DP.W[:,ind_set] = Sol[1:DP.n_full]
        for i_full=1:DP.n_full
            if DP.W[i_full,ind_set] != 0
                DP.Ws[i_full,ind_set] = Sym("W"*string(i_full)*"|"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"|"*string(ind_set))=>DP.W[i_full,ind_set])]
            end
        end
        DP.f[1:DP.n_aut,ind_set] = Sol[DP.n_full+1:DP.n_full+DP.n_aut]
        for i_aut=1:DP.n_aut
            if DP.f[i_aut,ind_set] != 0
                DP.fs[i_aut,ind_set] = Sym("f"*string(i_aut)*"|"*string(ind_set))
                DP.subs = [DP.subs;Dict(Sym("f"*string(i_aut)*"|"*string(ind_set))=>DP.f[i_aut,ind_set])]
            end
        end
    end

    """
    Auxiliary function for latex printing. Creates non commutative monoms so that they stay in
    the correct order when printing.
    This function is primarely intended for use inside the module, and not for exportation.
    """
    function non_commutative_monoms_for_latex_expressions(monom_number, monoms_exponents, vars)
        monom_string = ""
        for i = 1:length(vars)
            exp = monoms_exponents[monom_number][i]
            if exp != 0
                if exp == 1
                    monom_string *= "z_$(i)*"
                else
                    monom_string *= "z_$(i)^$(exp)*"
                end
            end
        end
        return sympy.Symbol(monom_string[1:end-1], commutative=False)
    end

    """
    Checks whether all of the terms in a sum have negative coefficients.
    This function is primarely intended for use inside the module, and not for exportation.
    """
    function all_negative_coeffs(expression)
        if expression.func == sympy.core.add.Add
            for arg in expression.args
                !sympy.core.function._coeff_isneg(arg) && return false
            end
        else
            !sympy.core.function._coeff_isneg(expression) && return false
        end
        return true
    end

    """
    Defines a consistent multiexponent ordering. For multiexponents m1 and m2,
    if deg1 < deg2, then m1 < m2. If deg1 = deg2, then the comparison is made
    looking at the degrees of each z_i in order. For example:
        - z1^3 * z2^4 * z3^1 > z1^3 * z2^3 * z3^3, because the z2 exponent is
        greater for the first monomial.
    """
    function multiexponent_lt(x, y)
        sum(x) < sum(y) && return true
        sum(y) < sum(x) && return false

        for i in eachindex(x)
            x[i] > y[i] && return true    
            y[i] > x[i] && return false
        end
        return false
    end

    """
    Auxiliary function to transform a given polynomial expression in its latex code. Input arguments are:
        - vars: A symbolic vector containing the variables of the polynomial expression.
        - expr_coeff: The polynomial coefficients.
        - expr_monoms: A vector of tuples with each monom's exponents.
        - latex_output: The string to which the latex code is appended.
    This function is primarely intended for use inside the module, and not for exportation.
    """
    function latex_code_for_polynomial_expression(vars::Vector{Sym}, expr_coeffs::Vector{Sym}, expr_monoms, latex_output::String)
        ordering = sortperm(expr_monoms, lt=multiexponent_lt)
        expr_coeffs = expr_coeffs[ordering]
        expr_monoms = expr_monoms[ordering]

        for j in eachindex(expr_coeffs)
            sign_flag = false
            monom = non_commutative_monoms_for_latex_expressions(j, expr_monoms, vars)
            # This sign flag is used as a way to make the i stay outside of fraction numerators
            if all_negative_coeffs(expr_coeffs[j])
                expr_coeffs[j] = -expr_coeffs[j]
                sign_flag = true
            end
            if !(expr_coeffs[j]/im).has(Sym(im))
                expr_coeffs[j] = sympy.Mul(Sym(im), expr_coeffs[j]/im, evaluate=False)
            end
            # In the next lines the Sym("a") is necessary so that -i is not displayed as \left( -i \right)
            if expr_coeffs[j] == 1
                latex_result = latexify(sympy.Add(Sym("a"), monom, evaluate=False), cdot = false)[3:end-1]
            else
                latex_result = latexify(sympy.Add(Sym("a"), sympy.Mul(expr_coeffs[j], monom, evaluate=False), evaluate=False), cdot = false)[3:end-1]
            end
            if sign_flag
                latex_result = " -"*latex_result[3:end]
            end
            if j == 1 && !sign_flag
                latex_result = latex_result[3:end]
            end
            latex_output *= latex_result
        end
        return latex_output
    end

    """
    Function to output the reduced dynamics on latex format. Input arguments are:\\
        - DP: The parametrisation data structure.\\
        - aexp: A multiexponent_struct defining the numbering of the monoms in the parametrisation.\\
        - substitutions: An array of dictionaries (or of arrays of dictionaries) containing substitutions
                         to be made on the symbolic variables.\\
        - output_file: A file to output the results. If it is not passed, results are displayed
                       on the terminal.\\
        - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                     If "w" overwrites existing files with the same name.
    """
    function reduced_dynamics_latex_output(DP::parametrisation_struct, aexp::multiexponent_struct, substitutions, output_file = nothing, file_mode::String = "a")
        poly_from_expr = sympy.polys.polytools.poly_from_expr

        println("Printing reduced dynamics")

        zₜ = sympy.zeros(DP.n_aut,1)[:,1]
        z = [Sym("z_$(i)") for i=1:DP.n_aut]
        for i_ord=1:length(DP.f[1,:])
            monom = prod(z .^ aexp.mat[:,i_ord])
            for i_var=1:DP.n_aut
                substituted = mysub([DP.f[i_var,i_ord]],DP.subs[end:-1:1])
                for substitute in substitutions
                    substituted = mysub(substituted, substitute)
                end
                zₜ[i_var:i_var] += substituted*monom
            end
        end
        
        latex_output = "\\begin{align}"
        for i in eachindex(zₜ)
            expr = poly_from_expr(zₜ[i], gens = z)
            expr_monoms = expr[1].monoms()
            expr_coeffs = expr[1].coeffs()
            latex_output *= "\n\\dot{z}_{$(i)} &="
            latex_output = latex_code_for_polynomial_expression(z, expr_coeffs, expr_monoms, latex_output)
            latex_output *= "\\\\"
        end
        latex_output = replace(latex_output, "I" => "\\mathit{i}")
        latex_output = latex_output[1:end-2] * "\n\\end{align}\n"
        
        if output_file === nothing
            println("Reduced dynamics:")
            println(latex_output)
        else
            open(output_file, file_mode) do file
                write(file, "Reduced dynamics:\n")
                write(file, latex_output)
            end 
        end
        println("Reduced dynamics printed")
    end

    """
    Function to output the nonlinear mappings on latex format. Input arguments are:\\
        - DP: The parametrisation data structure.\\
        - aexp: A multiexponent_struct defining the numbering of the monoms in the parametrisation.\\
        - substitutions: An array of dictionaries (or of arrays of dictionaries) containing substitutions
                         to be made on the symbolic variables.\\
        - output_file: A file to output the results. If it is not passed, results are displayed
                       on the terminal.\\
        - file_mode: The opening mode of the file. If "a", the output is appended at the end.
                     If "w" overwrites existing files with the same name.
    """
    function nonlinear_mappings_latex_output(DP::parametrisation_struct, aexp::multiexponent_struct, substitutions, output_file = nothing, file_mode::String = "a")
        poly_from_expr = sympy.polys.polytools.poly_from_expr

        println("Printing nonlinear mappings")

        u = sympy.zeros(DP.n_aut,1)[:,1]
        z = [Sym("z_$(i)") for i=1:DP.n_aut]
        for i_ord=1:length(DP.W[1,:])
            monom = prod(z .^ aexp.mat[:,i_ord])
            for i_var=1:DP.n_aut
                substituted = mysub([DP.W[i_var,i_ord]],DP.subs[end:-1:1])
                for substitute in substitutions
                    substituted = mysub(substituted, substitute)
                end
                u[i_var:i_var] += substituted*monom
            end
        end
        
        latex_output = "\\begin{align}"
        for i in eachindex(u)
            expr = poly_from_expr(u[i], gens = z)
            expr_monoms = expr[1].monoms()
            expr_coeffs = expr[1].coeffs()
            latex_output *= "\ny_{$(i)} &="
            latex_output = latex_code_for_polynomial_expression(z, expr_coeffs, expr_monoms, latex_output)
            latex_output *= "\\\\"
        end
        latex_output = replace(latex_output, "I" => "\\mathit{i}")
        latex_output = latex_output[1:end-2] * "\n\\end{align}\n"
        
        if output_file === nothing
            println("Nonlinear mappings:")
            println(latex_output)
        else
            open(output_file, file_mode) do file
                write(file, "Nonlinear mappings:\n")
                write(file, latex_output)
            end 
        end
        println("Nonlinear mappings printed")
    end
end




