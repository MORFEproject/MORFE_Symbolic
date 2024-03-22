using SymPy
using LinearAlgebra
push!(LOAD_PATH,joinpath(pwd(),"src"))
using MORFE_Symbolic

output_path = "./test/Mass_2_springs_conservative_unforced"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                System of equations definition                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_osc = 2                   # Number of oscillator equations
n_aux = 2                   # Number of auxiliary variables for quadratic recast
n_full = 2*n_osc + n_aux    # Number of variables on the full recast system

M, K, C = define_second_order_matrices(n_osc, mass = "unitary", stiffness = "diagonal", damping = nothing)
C[2,2] = C[1,1]

function F_nonlin(U,R,n_aux)
    F = sympy.zeros(n_osc,1)[:,1]

    # Coefficients depending on the eigenfrequencies
    ω = create_pos_vec("ω", n_osc)
    g¹₁₁ = Sym(3)/Sym(2)*ω[1]^2
    g²₁₁ = ω[2]^2/Sym(2);  g¹₁₂ = 2*g²₁₁;
    g¹₂₂ = ω[1]^2/Sym(2);  g²₁₂ = 2*g¹₂₂;
    g²₂₂ = Sym(3)/Sym(2)*ω[2]^2
    h¹₁₁₁ = (ω[1]^2+ω[2]^2)/Sym(2)
    h²₁₁₁ = 0; h¹₁₁₂ = 3*h²₁₁₁;
    h¹₁₂₂ = (ω[1]^2+ω[2]^2)/Sym(2); h²₁₁₂ = h¹₁₂₂;
    h¹₂₂₂ = 0; h²₁₂₂ = 3*h¹₂₂₂;
    h²₂₂₂ = (ω[1]^2+ω[2]^2)/Sym(2)
    
    F[1] = - (g¹₁₁*U[1]*U[1] + g¹₁₂*U[1]*U[2] + g¹₂₂*U[2]*U[2]) - 
             (h¹₁₁₁*U[1]*R[1] + h¹₁₁₂*U[2]*R[1] + h¹₁₂₂*U[1]*R[2] + h¹₂₂₂*R[2]*U[2])
    F[2] = - (g²₁₁*U[1]*U[1] + g²₁₂*U[1]*U[2] + g²₂₂*U[2]*U[2]) -
             (h²₁₁₁*U[1]*R[1] + h²₁₁₂*U[2]*R[1] + h²₁₂₂*U[1]*R[2] + h²₂₂₂*R[2]*U[2])
    return F
end

# Define static RHS vector 
C0 = sympy.zeros(n_full,1)[:,1]
# Define time-dependent RHS vectors
C⁺ₑₓₜ = sympy.zeros(n_full,1)[:,1]
C⁻ₑₓₜ = sympy.zeros(n_full,1)[:,1]

sys = define_system(n_aux, M, K, C, F_nonlin, C0, C⁺ₑₓₜ, C⁻ₑₓₜ)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                Definition of parametrisation                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_aut = 2                       # Number of complex coordinates in the autonomous part
n_nonaut = 0                    # Number of complex coordinates in the nonautonomous part
n_rom = n_aut + n_nonaut        # Number complex coordinates in the parametrisation
o = 5                           # Order of the expansion

λ₀ = create_gen_vec("λ",n_rom)  # Create a generic λ₀ vector (n_rom×1) 
# The first n_aut entries of λ₀ represent the master eigenvalues that will be chosen later
# remember to keep the same order here and there.Which means that if here two consecutive 
# entries represent a pair, also later they must the last n_nonaut entries of λ₀ represent ±imΩ.

conditions = [λ₀[2] =>-λ₀[1]]   # Resonance conditions
style = "CNF"                   # Parametrisation style
output_path *= "_" * style * "_output.txt"

DP, aexp = define_parametrisation(λ₀,n_rom,n_full,n_aut,n_osc,o,style,conditions)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                Computation of parametrisation                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

compute_order_0_parametrisation!(DP, aexp, sys)

eigenvalues_order = n_aux+2:-1:n_aux+1       # Ordering of the master eigenvalues
compute_order_1_parametrisation!(DP, aexp, sys, eigenvalues_order)

for p = 2:o
    compute_order_p_parametrisation!(DP, aexp, sys, p)
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                      Mathematica output                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# In this example Julia's simplification capabilities already give good outputs, the reason why 
# the next sections remain. However, to get the most simplified expressions, Mathematica is needed.
Mathematica_output(DP, aexp, output_path[1:end-11], "Output_Mathematica",
                    print_reduced_dynamics = true, print_nonlinear_mappings = true)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                         Substitutions                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
substitutions = []
substitutions!(DP, substitutions)
reduced_dynamics_substitutions!(DP, substitutions)
nonlinear_mappings_substitutions!(DP, substitutions)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                         Realification                                          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
omega, xi = backbone_CNF(DP, aexp)
amplitude = physical_amplitudes_CNF(DP, aexp)
cartesian_realification!(DP, aexp)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                           Printing                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
reduced_dynamics_latex_output(DP, aexp, output_path)
nonlinear_mappings_latex_output(DP, aexp, output_path)
backbone_output(omega, output_path)
physical_amplitudes_output(amplitude, output_path)