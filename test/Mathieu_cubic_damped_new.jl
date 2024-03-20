using SymPy
using LinearAlgebra
push!(LOAD_PATH,joinpath(pwd(),"src"))
using MORFE_Symbolic

# NOTE: This example is outside the scope of the function define_system. Because of that, the definition
# of the system of equations is done manually.

output_path = "./test/Mathieu_cubic_damped"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                System of equations definition                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

n_osc = 1                   # Number of oscillator equations
n_aux = 2                   # Number of auxiliary variables for quadratic recast
n_full = 2*n_osc + n_aux    # Number of variables on the full recast system

M, K, C = define_second_order_matrices(n_osc, mass = "unitary", stiffness = "diagonal", damping = "diagonal")

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
                                                    # and to define the parametric excitation
    return F
end

function RHS_Quad(Y)
    F = sympy.zeros(n_full,1)[:,1];
    U = Y[1:n_osc]                                  # First n_osc positions are the elements of U
    R = Y[2*n_osc+1:2*n_osc+n_aux]                  # Last n_aux positions are the auxiliary variables
    h = symbols("h", real = true)
    F[n_osc+1] = - h*U[1]*R[1] + U[1]*R[2]          # Nonlinear internal forces in quadratic form
    F[2*n_osc+1:2*n_osc+n_aux-1] = -U.^2            # Equations for the quadratic recast
    F[2*n_osc+n_aux] = 0                            # Last equation for the definition of the parametric excitation
    return F
end

# Define static RHS vector 
C0 = sympy.zeros(n_full,1)[:,1]
# Define time-dependent RHS vectors
C⁺ₑₓₜ = sympy.zeros(n_full,1)[:,1]
C⁻ₑₓₜ = sympy.zeros(n_full,1)[:,1]
# To have a κ*cosin(Ωt) excitation as R₂ for the definition of the parametric excitation
C⁺ₑₓₜ[2*n_osc+n_aux] = -1/Sym(2)*symbols("κ",positive=true)
C⁻ₑₓₜ[2*n_osc+n_aux] = -1/Sym(2)*symbols("κ",positive=true)

sys = system_struct(extract_Lin(RHS_Lin, n_full),
                    extract_Lin(LHS_Lin, n_full),
                    extract_Quad(RHS_Quad, n_full),
                    C0, C⁺ₑₓₜ, C⁻ₑₓₜ)

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

conditions = [λ₀[2] =>-λ₀[1],λ₀[4] =>-λ₀[3],λ₀[3] =>2*λ₀[1]]   # Resonance conditions
style = "CNF"                                                  # Parametrisation style
output_path *= "_" * style * "_output.txt"

DP, aexp = define_parametrisation(λ₀,n_rom,n_full,n_aut,n_osc,o,style,conditions)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                Computation of parametrisation                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

compute_order_0_parametrisation!(DP, aexp, sys)

eigenvalues_order = n_aux+2:-1:n_aux+1       # Ordering of the master eigenvalues
# Column vector containing the symbolic expressions of the nonautonomous eigenvalues
ω₁ = symbols("ω₁",positive = true)
nonautonomous_eigenvalues = [2*im*ω₁;-2*im*ω₁]
compute_order_1_parametrisation!(DP, aexp, sys, eigenvalues_order, nonautonomous_eigenvalues)

for p = 2:o
    compute_order_p_parametrisation!(DP, aexp, sys, p)
end

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                      Mathematica output                                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# In this example Julia's simplification capabilities can't get to the most compact form of the 
# equations. In such way, the substitutions are done in Mathematica.
Mathematica_output(DP, aexp, output_path[1:end-11], "Output_Mathematica",
                    print_reduced_dynamics = true, print_nonlinear_mappings = true)