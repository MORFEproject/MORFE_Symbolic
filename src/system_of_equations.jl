"""
Data structure defining the departing (full) system of equations:
    B.Yₜ = A.Y + (Q.Y).Y + C0 + C⁺ₑₓₜ exp(+ im Ω t) + C⁻ₑₓₜ exp(-im Ω t)
"""
mutable struct system_struct
    A::Matrix{Sym}
    B::Matrix{Sym}
    Q::Vector{Vector{Sym}}
    C0::Vector{Sym}
    C⁺ₑₓₜ::Vector{Sym}
    C⁻ₑₓₜ::Vector{Sym}
end

"""
Extracts the tensor Q from the quadratic function defining the system.
"""
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

"""
Extracts a matrix equivalent to the linear function LinFun.
"""
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

"""
Applies the Q tensor on vectors Y1 and Y2 and gives back the resulting vector.
"""
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

"""
Creates mass, stiffness and damping matrices relative to the second order problem.
Input parameters are:\\
    - n_osc: The size of the second order system (number of oscillator equations);\\
    - mass: A string defining the type of the mass matrix. Can be either unitary, diagonal or full;\\
    - stiffness: A string defining the type of the stiffness matrix. Can be either diagonal or full;\\
    - damping: A variable defining the type of the damping matrix. Can be either nothing (undamped system), diagonal or full.\\
An unitary mass matrix means that it is diagonal and with masses equal to one. When this option is chosen,
if the stiffness and damping matrices are diagonal, they are going to be consitituted by the frequencies ω^2 and by 2*ξ*ω,
respectively, with ξ the damping ratio.\\
Returns the matrices M, K and C.
"""
function define_second_order_matrices(n_osc; mass = "unitary", stiffness = "diagonal", damping = nothing)
    if mass == "unitary"
        M = diagm(sympy.ones(n_osc,1)[:,1])
    elseif mass == "diagonal"
        M = diagm(create_pos_vec("m",n_osc))
    elseif mass == "full"
        M = create_pos_matrix("m",n_osc)
    else 
        throw(ArgumentError("Mass matrix should be either unitary, diagonal or full!"))
    end

    if stiffness == "diagonal" && mass == "unitary"
        ω = create_pos_vec("ω",n_osc)
        K = diagm(ω.^2)
    elseif stiffness == "diagonal"
        K = diagm(create_real_vec("k",n_osc))
    elseif stiffness == "full"
        K = create_real_matrix("k",n_osc)
    else 
        throw(ArgumentError("Stiffness matrix should be either diagonal or full!"))
    end

    if damping === nothing
        C = diagm(sympy.zeros(n_osc,1)[:,1])
    else
        if damping == "diagonal" && mass == "unitary"
            ξ = create_pos_vec("ξ",n_osc)
            C = diagm(2*ξ.*ω)
        elseif damping == "diagonal"
            C = diagm(create_pos_vec("c",n_osc))
        elseif damping == "full"
            C = create_real_matrix("c",n_osc)
        else 
            throw(ArgumentError("Damping matrix should be either diagonal or full!"))
        end
    end

    return M, K, C
end

"""
Defines the system of equations structure. Input arguments are:\\
    - n_aux: The number of auxiliary variables for the quadratic recast;\\
    - M: The mass matrix;\\
    - K: The stiffness matrix;\\
    - C: The damping matrix;\\
    - F_nonlin: A function defining the nonlinear internal forces. Receives as input arguments the
    displacements U, the quadratic recast variables R and n_aux. Returns an array of length n_osc.\\
    - C0: An array defining the static forces;\\
    - C⁺ₑₓₜ and C⁻ₑₓₜ: Arrays defining the amplitudes of a single harmonic excitation of frequency Ω.
    The excitation is defined using complex variables, such that the + vector is associated with
    exp(iΩt) and the - one with exp(-iΩt).\\
    Returns the system of equations structure.
"""
function define_system(n_aux, M, K, C, F_nonlin, C0, C⁺ₑₓₜ, C⁻ₑₓₜ)
    n_osc = size(M)[1]
    n_full = 2*n_osc + n_aux
    
    function LHS_Lin(Yₜ)
        F = sympy.zeros(n_full,1)[:,1];
        Uₜ = Yₜ[1:n_osc]                                # First n_osc positions are the elements of U
        Vₜ = Yₜ[n_osc+1:2*n_osc]                        # Following n_osc positions are the elements of V
        F[1:n_osc] = M*Uₜ                               # First n_osc equations is the LHS of M*Uₜ = M*V
        F[n_osc+1:2*n_osc] = M*Vₜ                       # Following n_osc equations is the LHS of M*Vₜ = -C*V -K*U 
        # Last n_aux equations are for the quadratic recast so their LHS is zero
        return F
    end

    function RHS_Lin(Y)
        F = sympy.zeros(n_full,1)[:,1];
        U = Y[1:n_osc]                                  # First n_osc positions are the elements of U
        V = Y[n_osc+1:2*n_osc]                          # Following n_osc positions are the elements of V
        R = Y[2*n_osc+1:2*n_osc+n_aux]                  # Last n_aux positions are the auxiliary variables
        F[1:n_osc] = M*V                                # First n_osc equations is the RHS of M*Uₜ = M*V
        F[n_osc+1:2*n_osc] = -C*V-K*U                   # Following n_osc equations is the RHS of M*Vₜ = -C*V -K*U 
        F[2*n_osc+1:2*n_osc+n_aux] = R                  # Last n_aux equations are for the quadratic recast
        return F
    end

    function RHS_Quad(Y)
        F = sympy.zeros(n_full,1)[:,1];
        U = Y[1:n_osc]                                  # First n_osc positions are the elements of U
        R = Y[2*n_osc+1:2*n_osc+n_aux]                  # Last n_aux positions are the auxiliary variables
        F[n_osc+1:2*n_osc] = F_nonlin(U,R,n_aux)        # Nonlinear internal forces in quadratic form
        F[2*n_osc+1:2*n_osc+n_aux] = -U.^2              # Last n_aux equations are for the quadratic recast
        return F
    end

    sys = system_struct(extract_Lin(RHS_Lin, n_full),
                        extract_Lin(LHS_Lin, n_full),
                        extract_Quad(RHS_Quad, n_full),
                        C0, C⁺ₑₓₜ, C⁻ₑₓₜ)

    return sys
end