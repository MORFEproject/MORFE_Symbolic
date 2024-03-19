using SymPy
using LinearAlgebra
push!(LOAD_PATH,joinpath(pwd(),"src"))
using MORFE_Symbolic
using Combinatorics

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                System of equations definition                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_osc = 1                   # Number of oscillator equations
n_aux = 1                   # Number of auxiliary variables for quadratic recast
n_full = 2*n_osc + n_aux    # Number of variables on the full recast system

M, K, C = define_second_order_matrices(n_osc, mass = "unitary", stiffness = "diagonal")

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
n_nonaut = 0                    # Number of complex coordinates in the nonautonomous part
n_rom = n_aut + n_nonaut        # Number complex coordinates in the parametrisation
o = 11                          # Order of the expansion

λ₀ = create_gen_vec("λ",n_rom)  # Create a generic λ₀ vector (n_rom×1) 
# The first n_aut entries of λ₀ represent the master eigenvalues that will be chosen later
# remember to keep the same order here and there.Which means that if here two consecutive 
# entries represent a pair, also later they must the last n_nonaut entries of λ₀ represent ±imΩ.

conditions = [λ₀[2] =>-λ₀[1]]   # Resonance conditions
style = "CNF"                   # Parametrisation style
DP, aexp = define_parametrisation(λ₀,n_rom,n_full,n_aut,n_osc,o,style,conditions)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                Computation of parametrisation                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

compute_order_0_parametrisation!(DP, aexp, sys)

eigenvalues_order = n_aux+2:-1:n_aux+1       # Ordering of the master eigenvalues
# Column vector containing the symbolic expressions of the nonautonomous eigenvalues
nonautonomous_eigenvalues = [im*symbols("Ω",positive=true);-im*symbols("Ω",positive=true)]
compute_order_1_parametrisation!(DP, aexp, sys, eigenvalues_order, nonautonomous_eigenvalues)

for p = 2:o
    compute_order_p_parametrisation!(DP, aexp, sys, p)
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                         Substitutions                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
ξ = create_pos_vec("ξ",n_osc)
δ = create_pos_vec("δ",n_osc)
substitutions = [[Dict(sqrt(ξ[i]^2 - 1)=>im*δ[i]) for i=1:n_osc], [Dict(2*ξ[i]^3 - 2*ξ[i] =>-2*ξ[i]δ[i]^2) for i=1:n_osc]]
substitutions!(DP, substitutions)
reduced_dynamics_substitutions!(DP, substitutions)
nonlinear_mappings_substitutions!(DP, substitutions)
modal_coordinates_from_physical_coordinates!(DP, aexp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                         Realification                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
omega, xi = backbone_CNF(DP, aexp)
amplitude = physical_amplitudes_CNF(DP, aexp)
cartesian_realification!(DP, aexp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                           Printing                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
reduced_dynamics_latex_output(DP, aexp, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt")
nonlinear_mappings_latex_output(DP, aexp, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt")
nonlinear_mappings_latex_output(DP, aexp, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt", result = "modal")
reduced_dynamics_latex_output(DP, aexp, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt", normal_coordinate = 'a', result = "real")
nonlinear_mappings_latex_output(DP, aexp, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt", normal_coordinate = 'a', result = "real")
backbone_output(omega, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt")
nonlinear_damping_output(omega, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt")
physical_amplitudes_output(amplitude, "./test/Duffing_cubic_conservative_unforced_CNF_output.txt")