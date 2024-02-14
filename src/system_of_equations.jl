"""
Data structure defining the departing (full) system of equations:
    A.Yₜ = B.Y + (Q.Y).Y + C0 + C⁺ₑₓₜ exp(+ im Ω t) + C⁻ₑₓₜ exp(-im Ω t)
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

function define_system(n_osc, n_aux, conservative)
    M = diagm(sympy.ones(n_osc,1)[:,1])
    ω = create_pos_vec("ω",n_osc)
    K = diagm(ω.^2)

    ξ = conservative ? sympy.zeros(n_osc,1)[:,1] : create_pos_vec("ξ",n_osc)
    ζ = 2*ξ.*ω
    C = diagm(ζ)

end