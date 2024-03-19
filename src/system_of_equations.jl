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

function define_nonlinar_tensors(n_osc; quadratic = nothing, cubic = nothing)
    if quadratic == "full"
        G = create_real_3_tensor("g",n_osc)
    elseif quadratic == "diagonal"
        G = create_real_vec("g",n_osc)
    elseif quadratic === nothing
        G = nothing
    else 
        throw(ArgumentError("Quadratic nonlinear tensor should be either full, diagonal or nothing!"))
    end

    if cubic == "full"
        H = create_real_4_tensor("h",n_osc)
    elseif cubic == "diagonal"
        H = create_real_vec("h",n_osc)
    elseif cubic === nothing
        H = nothing
    else 
        throw(ArgumentError("Cubic nonlinear tensor should be either full, diagonal or nothing!"))
    end

    return G, H
end

# function quadratic_stiffness(n,G,Y1,Y2,symmetry)
#     if G === nothing
#         return sympy.zeros(n,1)[:,1]
#     end
#     if typeof(G) == Vector{Sym}
#         return G.*Y1.*Y2
#     end

#     res = sympy.zeros(n,1)[:,1]
#     visited = Set()
#     for i = 1:n
#         for j = 1:n
#             for k = 1:n
#                 if symmetry && !([i,j,k] in visited)
#                     val = G[i,j,k]
#                     for index in permutations([i,j,k])
#                         push!(visited, index)
#                         G[index[1],index[2],index[3]] = val
#                     end
#                 end
#                 res[i] += G[i,j,k]*Y1[j]*Y2[k]
#             end
#         end
#     end

#     return res 
# end

# function cubic_stiffness(n,H,Y1,Y2,Y3,symmetry)
#     if H === nothing
#         return sympy.zeros(n,1)[:,1]
#     end
#     if typeof(H) == Vector{Sym}
#         return H.*Y1.*Y2.*Y3
#     end

#     res = sympy.zeros(n,1)[:,1]
#     visited = Set()
#     for i = 1:n
#         for j = 1:n
#             for k = 1:n
#                 for l = 1:n
#                     if symmetry && !([i,j,k,l] in visited)
#                         val = G[i,j,k,l]
#                         for index in permutations([i,j,k,l])
#                             push!(visited, index)
#                             H[index[1],index[2],index[3],index[4]] = val
#                         end
#                     end
#                     res[i] += H[i,j,k]*Y1[j]*Y2[k]*Y3[l]
#                 end
#             end
#         end
#     end

#     return res 
# end

function define_system(n_aux, M, K, C, F_nonlin, C0, C⁺ₑₓₜ, C⁻ₑₓₜ)
    n_osc = size(M)[1]
    n_full = 2*n_osc + n_aux
    
    function LHS_Lin(Yₜ)
        F = sympy.zeros(n_full,1)[:,1];
        Uₜ = Yₜ[1:n_osc]                    # first n_osc positions is U
        Vₜ = Yₜ[n_osc+1:2*n_osc]      # second n_osc positions is V
        # first n_osc equations is M*Uₜ = M*V
        F[1:n_osc] = M*Uₜ
        # second n_osc equations is M*Vₜ = -C*V -K*U ...
        F[n_osc+1:2*n_osc] = M*Vₜ
        # last n_aux equations are the algebraic ones defining the auxiliary variables
        # so they do not have a LHS
        return F
    end

    function RHS_Lin(Y)
        F = sympy.zeros(n_full,1)[:,1];
        U = Y[1:n_osc]                                    # first n_osc positions is U
        V = Y[n_osc+1:2*n_osc]                      # second n_osc positions is V
        R = Y[2*n_osc+1:2*n_osc+n_aux]       # last n_aux positions are the auxiliary variables
        # first n_osc equations is M*Uₜ = M*V
        F[1:n_osc] = M*V
        # second n_osc equations is M*Vₜ = -C*V -K*U ...
        F[n_osc+1:2*n_osc] = -C*V-K*U
        # last n_aux equations are the algebraic ones defining the auxiliary variables
        F[2*n_osc+1:2*n_osc+n_aux] = R
        return F
    end

    function RHS_Quad(Y)
        F = sympy.zeros(n_full,1)[:,1];
        U = Y[1:n_osc]                      # first n_osc positions is U
        R = Y[2*n_osc+1:2*n_osc+n_aux]      # last n_aux positions are the auxiliary variables
        # assign to the second n_osc equations
        F[n_osc+1:n_osc+n_aux] = F_nonlin(U,R,n_aux)
        # last n_aux equations are the algebraic ones defining the auxiliary variables
        F[2*n_osc+1:2*n_osc+n_aux] = -U.^2
        return F
    end

    sys = system_struct(extract_Lin(   RHS_Lin,    n_full),
                            extract_Lin(   LHS_Lin,    n_full),
                            extract_Quad(RHS_Quad,n_full),
                            C0,C⁺ₑₓₜ,C⁻ₑₓₜ)

    return sys
end