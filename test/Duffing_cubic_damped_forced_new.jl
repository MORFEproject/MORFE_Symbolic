using SymPy
using LinearAlgebra
push!(LOAD_PATH,joinpath(pwd(),"src"))
using MORFE_Symbolic

output_path = "./test/Duffing_cubic_damped_forced"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                System of equations definition                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_osc = 1                   # Number of oscillator equations
n_aux = 1                   # Number of auxiliary variables for quadratic recast
n_full = 2*n_osc + n_aux    # Number of variables on the full recast system

M, K, C = define_second_order_matrices(n_osc, mass = "unitary", stiffness = "diagonal", damping = "diagonal")

function F_nonlin(U,R,n_aux)
    F = sympy.zeros(n_osc,1)[:,1]
    h = symbols("h", real=true)
    F[1] = -h*U[1]*R[1]
    return F
end

# Define static RHS vector 
C0 = sympy.zeros(n_full,1)[:,1]
# Define time-dependent RHS vectors
C⁺ₑₓₜ = sympy.zeros(n_full,1)[:,1]
C⁻ₑₓₜ = sympy.zeros(n_full,1)[:,1]
# To have a κ*cosin(Ωt) excitation on U₁
C⁺ₑₓₜ[n_osc+1] = 1/Sym(2)*symbols("κ",positive=true)
C⁻ₑₓₜ[n_osc+1] = 1/Sym(2)*symbols("κ",positive=true)

sys = define_system(n_aux, M, K, C, F_nonlin, C0, C⁺ₑₓₜ, C⁻ₑₓₜ)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                Definition of parametrisation                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_aut = 2                       # Number of complex coordinates in the autonomous part
n_nonaut = 2                    # Number of complex coordinates in the nonautonomous part
n_rom = n_aut + n_nonaut        # Number complex coordinates in the parametrisation
o = 3                           # Order of the expansion

λ₀ = create_gen_vec("λ",n_rom)  # Create a generic λ₀ vector (n_rom×1) 
# The first n_aut entries of λ₀ represent the master eigenvalues that will be chosen later
# remember to keep the same order here and there.Which means that if here two consecutive 
# entries represent a pair, also later they must the last n_nonaut entries of λ₀ represent ±imΩ.

conditions = [λ₀[2] =>-λ₀[1],λ₀[4] =>-λ₀[3]]   # Resonance conditions
resonance = "Primary"                          # Either: Out-of-resonace, Primary, Superharmonic, Subharmonic
if resonance == "Primary"
    push!(conditions, λ₀[3] => λ₀[1])
elseif resonance == "Superharmonic"
    push!(conditions, 3*λ₀[3] => λ₀[1])
elseif resonance == "Subharmonic"
    push!(conditions, λ₀[3] => 3*λ₀[1])
end

style = "CNF"                                  # Parametrisation style
output_path *= "_" * style * "_output.txt"

DP, aexp = define_parametrisation(λ₀,n_rom,n_full,n_aut,n_osc,o,style,conditions)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                Computation of parametrisation                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

compute_order_0_parametrisation!(DP, aexp, sys)

eigenvalues_order = n_aux+2:-1:n_aux+1       # Ordering of the master eigenvalues
# Column vector containing the symbolic expressions of the nonautonomous eigenvalues
ω₁ = symbols("ω₁",positive = true)
if resonance == "Out-of-resonance"
    nonautonomous_eigenvalues = [im*symbols("Ω",positive=true);-im*symbols("Ω",positive=true)]
elseif resonance == "Primary"
    nonautonomous_eigenvalues = [im*ω₁;-im*ω₁]
elseif resonance == "Superharmonic"
    nonautonomous_eigenvalues = [im*ω₁/3;-im*ω₁/3]
elseif resonance == "Subharmonic"
    nonautonomous_eigenvalues = [3*im*ω₁;-3*im*ω₁]
end
compute_order_1_parametrisation!(DP, aexp, sys, eigenvalues_order, nonautonomous_eigenvalues)

for p = 2:o
    compute_order_p_parametrisation!(DP, aexp, sys, p)
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                      Mathematica output                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# In this example Julia's simplification capabilities can't get to the most compact form of the 
# equations. In such way, the substitutions are done in Mathematica.
Mathematica_output(DP, aexp, output_path[1:end-11]*"/"*resonance, "Output_Mathematica",
                    print_reduced_dynamics = true, print_nonlinear_mappings = true)