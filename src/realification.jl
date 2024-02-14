using SymPy

# At present limited to one master mode with one harmonic excitation 
function polar_realification(DP::parametrisation_struct,aexp::multiexponent_struct)
    get_re_im = sympy.core.expr.Expr.as_real_imag
    trigsimp = sympy.core.expr.Expr.trigsimp
    Wild = sympy.core.symbol.Wild
    collect = sympy.core.expr.Expr.collect

    t1 = time_ns()
    println("Polar realification started")

    θ = symbols("θ", real=true); ρ = symbols("ρ", real=true); Ω = symbols("Ω", real=true)
    real = 0; imaginary = 0

    z = [ρ/2*exp(im*θ), ρ/2*exp(-im*θ)]
    if DP.n_aut != DP.n_rom
        append!(z,[exp(im*Ω)/2, exp(-im*Ω)/2])
    end

    for i_ord in 1:length(DP.f[1,:])
        term = DP.f[1,i_ord]*prod(z .^ aexp.mat[:,i_ord])/exp(im*θ)
        term = get_re_im(term.expand(complex = true))
        real += simplify(trigsimp(term[1]))
        imaginary += simplify(trigsimp(term[2]))
    end
    m = Wild("m"); n = Wild("n")
    real_parts = collect(expand(2*real),(sin(Ω-θ),cos(Ω-θ)), evaluate = false)
    imaginary_parts = collect(expand(2*imaginary),(sin(Ω-θ),cos(Ω-θ)), evaluate = false)

    real = 0; imaginary = 0;
    for key in keys(real_parts)
        real = sympy.Add(real, simplify(real_parts[key])*key, evaluate=False)
    end
    for key in keys(imaginary_parts)
        imaginary = sympy.Add(imaginary, simplify(imaginary_parts[key])*key, evaluate=False)
    end
    real = sympy.apart(real, ρ)
    imaginary = sympy.apart(imaginary, ρ)

    t2 = time_ns()
    println("Polar realification ended")
    println("Elapsed time: $((t2-t1)/1.0e9) s")
    println("")

    return real, imaginary
end

function cartesian_realification!(DP::parametrisation_struct, aexp::multiexponent_struct, n_aux::Int64)
    t1 = time_ns()
    println("Cartesian realification started")

    for i_set in 1:DP.n_sets
        Av = aexp.mat[:,i_set]
        p = sum(Av)
        Iv = zeros(Int64,p)
        counter = 0
        for j in 1:DP.n_rom 
            Iv[counter+1:counter+Av[j]] .= j
            counter += Av[j]
        end
        pos = 1
        coeff = Sym(1)
        Av .= 0
        recursive_C2R!(Iv, p, pos, i_set, Av, coeff, DP, aexp)
    end

    nzhalf = Int(DP.n_aut/2)
    nzfhalf=Int((DP.n_rom-DP.n_aut)/2)
    for i in 1:nzhalf
        DP.fr[i+nzhalf,:] = 2*imag(DP.fr[i,:])
        DP.fr[i,:] = 2*real(DP.fr[i,:])
    end
    for i in DP.n_aut+1:DP.n_aut+nzfhalf
        DP.fr[i+nzfhalf,:] = imag(DP.fr[i,:])
        DP.fr[i,:] = real(DP.fr[i,:])
    end
    DP.Wr = real(DP.Wr)

    t2 = time_ns()
    println("Cartesian realification ended")
    println("Elapsed time: $((t2-t1)/1.0e9) s")
    println("")
end

function recursive_C2R!(Iv::Vector{Int64}, p::Int64, pos::Int64, posinit::Int64, Av::Vector{Int64},
    coeff::Sym, DP::parametrisation_struct, aexp::multiexponent_struct)

    nzhalf=Int(DP.n_aut/2)
    nzfhalf=Int((DP.n_rom-DP.n_aut)/2)                           
    if pos == (p+1)
        pos1 = aexp.get(Av)
        DP.Wr[:,pos1] += coeff*DP.W[:,posinit]
        DP.fr[:,pos1] += coeff*DP.f[:,posinit]
    else
        Av1 = Av[:]   # new vectors
        Av2 = Av[:]
        iz = Iv[pos]    # z var 

        if iz <= nzhalf
            coeff1 = coeff/2  
            Av1[iz] += 1
            coeff2 = im*coeff/2
            Av2[iz+nzhalf] += 1
        elseif iz <= DP.n_aut
            coeff1 = coeff/2
            Av1[iz-nzhalf] += 1
            coeff2 = -im*coeff/2
            Av2[iz] += 1
        elseif iz <= DP.n_aut+nzfhalf  # to have cos and sin the 1/2 must avoided
            coeff1 = coeff
            Av1[iz] += 1
            coeff2 = im*coeff
            Av2[iz+nzfhalf] += 1
        else    
            coeff1 = coeff
            Av1[iz-nzfhalf] += 1
            coeff2 = -im*coeff
            Av2[iz] += 1
        end     

        pos += 1
        recursive_C2R!(Iv, p, pos, posinit, Av1, coeff1, DP, aexp)
        recursive_C2R!(Iv, p, pos, posinit, Av2, coeff2, DP, aexp)
    end

end

# function recursive_C2R!(Iv::Vector{Int64}, p::Int64, pos::Int64, posinit::Int64, Av::Vector{Int64},
#     coeff::ComplexF64, DP::parametrisation_struct, aexp::multiexponent_struct)

#     nzhalf=Int(DP.n_aut/2)
#     nzfhalf=Int((DP.n_rom-DP.n_aut)/2)                           
#     if pos == (p+1)
#         pos1 = aexp.get(Av)
#         DP.Wr[:,pos1] += coeff*DP.W[:,posinit]
#         DP.fr[:,pos1] += coeff*DP.f[:,posinit]
#     else
#         Av1 = Av[:]   # new vectors
#         Av2 = Av[:]
#         iz = Iv[pos]    # z var 

#         if iz <= nzhalf
#             coeff1 = -0.5*coeff  
#             Av1[iz] += 1
#             coeff2 = 0.5*im*coeff
#             Av2[iz+nzhalf] += 1
#         elseif iz <= DP.n_aut
#             coeff1 = -0.5*coeff
#             Av1[iz-nzhalf] += 1
#             coeff2 = -0.5*im*coeff
#             Av2[iz] += 1
#         elseif iz <= DP.n_aut+nzfhalf  # to have cos and sin the 1/2 must avoided
#             coeff1 = -coeff
#             Av1[iz] += 1
#             coeff2 = -im*coeff
#             Av2[iz+nzfhalf] += 1
#         else    
#             coeff1 = coeff
#             Av1[iz-nzfhalf] += 1
#             coeff2 = -im*coeff
#             Av2[iz] += 1
#         end     

#         pos += 1
#         recursive_C2R!(Iv, p, pos, posinit, Av1, coeff1, DP, aexp)
#         recursive_C2R!(Iv, p, pos, posinit, Av2, coeff2, DP, aexp)
#     end

# end

function backbone_CNF(DP::parametrisation_struct, aexp::multiexponent_struct, omega = nothing)
    get_re_im = sympy.core.expr.Expr.as_real_imag
    trigsimp = sympy.core.expr.Expr.trigsimp

    t1 = time_ns()
    println("Backbone calculation started")

    θ = symbols("θ", real=true); ρ = symbols("ρ", real=true)

    z = [ρ/2*exp(im*θ), ρ/2*exp(-im*θ)]
    omega_rho = sympy.zeros(1,DP.order)
    xi_rho = sympy.zeros(1,DP.order)
    for i_set in 2:length(DP.f[1,:]) # Disregarding constant terms
        order = sum(aexp.mat[:,i_set])
        term = 2*DP.f[1,i_set]*prod(z .^ aexp.mat[:,i_set])/exp(im*θ)/ρ
        term = get_re_im(term.expand(complex = true))
        real = simplify(trigsimp(term[1]))
        imag = simplify(trigsimp(term[2]))
        if real != 0
            println("Caution! There seems to be non-imaginary coefficients on the backbone calculation.")
        end
        omega_rho[order] += imag
        xi_rho[order] -= real
    end

    if omega === nothing
        omega = omega_rho[1]
    end
    xi_rho = xi_rho./omega

    t2 = time_ns()
    println("Backbone calculation ended")
    println("Elapsed time: $((t2-t1)/1.0e9) s")
    println("")
    return omega_rho, xi_rho
end

# Function only for one master mode without damping
function physical_amplitudes_CNF(DP::parametrisation_struct, aexp::multiexponent_struct)
    get_re_im = sympy.core.expr.Expr.as_real_imag
    trigsimp = sympy.core.expr.Expr.trigsimp

    t1 = time_ns()
    println("Physical amplitudes calculation started")

    θ = symbols("θ", real=true); ρ = symbols("ρ", positive=true)

    z = [ρ/2*exp(im*θ), ρ/2*exp(-im*θ)]
    real_part = sympy.zeros(1,DP.order); imaginary_part = sympy.zeros(1,DP.order)
    for i_set in 2:length(DP.W[1,:]) # Disregarding constant terms
        order = sum(aexp.mat[:,i_set])
        term = DP.W[1,i_set]*ρ^order/2^order
        term = get_re_im(term.expand(complex = true))
        real_part[order] += term[1]
        if aexp.mat[1,i_set] >= aexp.mat[2,i_set]
            imaginary_part[order] += term[2]
        else
            imaginary_part[order] -= term[2]
        end
    end

    t2 = time_ns()
    println("Physical amplitudes calculation ended")
    println("Elapsed time: $((t2-t1)/1.0e9) s")
    println("")

    if sum(real_part) == 0
        for order in 1:DP.order
            imaginary_part[order] = simplify(imaginary_part[order])
        end
        return imaginary_part
    elseif sum(imaginary_part) == 0
        for order in 1:DP.order
            real_part[order] = simplify(real_part[order])
        end
        return real_part
    else
        throw(ErrorException("There is an imaginary part on the displacement field."))
    end
end

function modal_coordinates_from_physical_coordinates!(DP::parametrisation_struct, eigenvecs, n_osc)
    red_eigenvecs = eigenvecs[1:2*n_osc,:]
    for i in 1:length(DP.W[1,:])
        DP.Wmodal[1:2*n_osc,i] = red_eigenvecs \ DP.W[1:2*n_osc,i]
    end
end