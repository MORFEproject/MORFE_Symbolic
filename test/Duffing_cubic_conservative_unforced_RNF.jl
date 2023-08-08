using SymPy
using LinearAlgebra
push!(LOAD_PATH,joinpath(pwd(),"src"))
using MORFE_Symbolic


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                            Definition of original system                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#= 
Example of system definition with quadratic recast and first order

M Uₜ = M V
M Vₜ = - K U - C V - H₁₁₁ R₁ U₁  - F0 - F⁺ₑₓₜ exp(+ im Ω t) + F⁻ₑₓₜ exp(-im Ω t)
    0 = R₁ - U₁²


Y[1] = U₁
Y[2] = V₁
Y[3] = R₁

The user can write the equations in the functions:
LHS_Lin(Yₜ) = RHS_Lin(Y) + RHS_Quad(Y) + C0 + C⁺ₑₓₜ exp(+ im Ω t) + C⁻ₑₓₜ exp(-im Ω t)

Upon extraction, the system will read:
A.Yₜ = B.Y + (Q.Y).Y + C0 + C⁺ₑₓₜ exp(+ im Ω t) + C⁻ₑₓₜ exp(-im Ω t)

Q(Y)
(Q.Y1).Y2
Q(Y1,Y2)

Where the matrices A and B, the tensor Q, and the vectors C0,C⁺ₑₓₜ,C⁻ₑₓₜ
    will be stored in sys as: 
sys.A::Matrix{Sym}
sys.B::Matrix{Sym}
sys.Q::Dict{Any,Any}
sys.C0::Vector{Sym}
sys.C⁺ₑₓₜ::Vector{Sym}
sys.C⁻ₑₓₜ::Vector{Sym}

 =#


# input any mass matrix:
#
# m = symbols("m",positive = true)
# M = m
#
# or simply define an identity mass matrix of size n_osc:
#
n_osc = 1 # size of the original system in oscillatory form
M = diagm(sympy.ones(n_osc,1)[:,1])
#
# define a generic stiffness matrix:
# k = symbols("k",positive = true)
# K = k
#
# or simply define a diagonal matrix with entries ωⱼ^2:
# n_osc = size(M)[1]
ω = create_pos_vec("ω",n_osc)
K = diagm(ω.^2)
#
# if nonconservative
# 
# create a generic damping matrix:
# c = symbols("c",positive = true)
# C = c
#
# or simply create a diagonalised damping matrix
# generic diagonal damping:
ξ = create_pos_vec("ξ",n_osc)
ζ = 2*ξ.*ω
C = diagm(ζ)
# for the sake of readability, it is useful to specify that each oscillator is underdamped
# which means that the quantity   δⱼ := √(1-ξⱼ^2) is positive
# definition of δ = √(1-ξ.^2) will be used later for simplification
δ = create_pos_vec("δ",n_osc)

# the total size of the DAE system will be 
# the size of the oscillatory system in first order form (2*n_osc)
# plus the number of algebraic equations needed for the quadratic recast
# is n_osc = 1 and the nonlinearity is cubic, 
# only one auxiliary variable is needed (R₁ = U₁^2)
n_aux = 1
n_full = 2*n_osc+n_aux

# define the LHS as a function LHS_Lin(Yₜ)
# the matrix A such that LHS_Lin(Yₜ) = A.Yₜ
# will be automatically extracted later
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

# define the Linear RHS as a function RHS_Lin(Y)
# the matrix B such that RHS_Lin(Y) = B.Y
# will be automatically extracted later
function RHS_Lin(Y)
    F = sympy.zeros(n_full,1)[:,1];
    U = Y[1:n_osc]                                    # first n_osc positions is U
    V = Y[n_osc+1:2*n_osc]                      # second n_osc positions is V
    R = Y[2*n_osc+1:2*n_osc+n_aux]       # last n_aux positions are the auxiliary variables
    # first n_osc equations is M*Uₜ = M*V
    F[1:n_osc] = M*V
    # second n_osc equations is M*Vₜ = -C*V -K*U ...
    #F[n_osc+1:2*n_osc] = -C*V-K*U
    # NB : pour rendre le système conservatif, ici j'ai viré l'amortissement.
    F[n_osc+1:2*n_osc] = -K*U
    # last n_aux equations are the algebraic ones defining the auxiliary variables
    F[2*n_osc+1:2*n_osc+n_aux] = R
    return F
end

# define the Quadratic RHS as a function RHS_Quad(Y)
# the sparse tensor Q such that RHS_Quad(Y) = (Q.Y).Y
# will be automatically extracted later
function RHS_Quad(Y)
    F = sympy.zeros(n_full,1)[:,1];
    U = Y[1:n_osc]                        # first n_osc positions is U
    R = Y[2*n_osc+1:2*n_osc+n_aux]      # last n_aux positions are the auxiliary variables
    # define only cubic nonlinearity
    h = Sym("h")
    # assign to the second n_osc equations
    # la ligne ci-dessous pour quad et cub: commentée
    #F[n_osc+1] = -  h*U[1]*R[1] - g*U[1]*U[1]
    #pour commencer on ne garde que le cubique
    F[n_osc+1] = -  h*U[1]*R[1]
    # last n_aux equations are the algebraic ones defining the auxiliary variables
    F[2*n_osc+1:2*n_osc+n_aux] = -U.^2
    return F
end

# define static RHS vector 
C0 = sympy.zeros(n_full,1)[:,1]

# define time-dependent RHS vectors
C⁺ₑₓₜ = sympy.zeros(n_full,1)[:,1]
C⁻ₑₓₜ = sympy.zeros(n_full,1)[:,1]
# to have a κ*cosin(Ωt) excitation on U₁
C⁺ₑₓₜ[n_osc+1] = -1/Sym(2)*symbols("κ",positive=true)
C⁻ₑₓₜ[n_osc+1] = -1/Sym(2)*symbols("κ",positive=true)
# to have a κ*sin(Ωt) excitation on U₁
# C⁺ₑₓₜ[n_osc+1] = + im/Sym(2)*symbols("κ",positive=true)
# C⁻ₑₓₜ[n_osc+1] = -  im/Sym(2)*symbols("κ",positive=true)

# the function extract_Lin extracts the A and B matrices
# from the user defined functions LHS_Lin and RHS_Lin
# the function extract_Quad extracts the sparse tensor Q
# from the user defined function RHS_Quad
# the structure sys contains: sys.A, sys.B, sys.C0, sys.Q, sys.C⁺ₑₓₜ, sys.C⁻ₑₓₜ
sys = system_struct(extract_Lin(   LHS_Lin,    n_full),
                            extract_Lin(   RHS_Lin,    n_full),
                            extract_Quad(RHS_Quad,n_full),
                            C0,C⁺ₑₓₜ,C⁻ₑₓₜ)




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                            Definition of parametrisation                                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# parametrisation variables 
# number of complex coordinates in the autonomous part
n_aut = 2
#
# number of complex coordinates in the nonautonomous part
# since the case of multiple forcing frequencies is not treated, the problem is
# either autonomous or nonautonomous with a single frequency and 2 complex coordinates
# for this reason, n_nonaut is either 0 or 2
n_nonaut = 0
#
# total number complex coordinates in the parametrisation
n_rom = n_aut + n_nonaut
#
# order of the expansion
o = 9
#
# initialise aexp
# this is a structure containing information about all the sets
# relative to each monomial in z, normal coordinates
# ATTENTION: the normal coordinates are both the ones relating
# to the autonomous part, both the ones for the nonautonomous
# the latter are automatically defined by:
# z[n_aut+1] = exp(+imΩt)
# z[n_aut+2] = exp(-imΩt)
aexp = init_multiexponent_struct(n_rom,o)
# aexp.mat is a (n_rom×n_sets) matrix
# aexp.n_sets is the total number of monomials
#
# each row of aexp.mat contains the exponents relating to a given set
# in a given row, each column k represent the exponent relating to the coordinate zₖ
# so that the I-th monomial will be equal to:
# ∑ₖ₌₁ⁿʳᵒᵐ   z[k]^aexp.mat[k,I]
# try typing in the terminal:
# prod(create_gen_vec("z",2).^aexp.mat,dims=1)
# to see how the monomials are organised

# initialise DP
# this is a structure that will contain the solution of each step
# of the parametrisation method
# here it is only initialised with zeros
DP = init_parametrisation_struct(n_full,n_rom,aexp.n_sets,n_aut)
# DP.W is a (n_full×n_sets) matrix whose colums contain the mapping 
# relating to each monomial 
# Y = ∑  DP.W[:,I]*z^aexp[I,:]
# DP.f is a (n_rom×n_sets) matrix whose colums contain the reduced dynamics
# ∂ₜ zₖ =  ∑  DP.f[k,I]*z^aexp[I,:]
# DP.σ is a (1×n_sets) vector containing the summation of λ relating to each monomial
# DP.σ[I] = λ*aexp[I,:]

# resonance conditions:
# create a generic λ₀ vector (n_rom×1) 
λ₀ = create_gen_vec("λ₀",n_rom)
# the first n_aut entries of λ₀ represent the master eigenvalues that will be chosen later
# remember to keep the same order here and there
# which means that if here two consecutive entries represent a pair, also later they must
# the last n_nonaut entries of λ₀ represent ±imΩ
#
# let us assume that the rom is an unforced oscillator 
# then the conditions of near resonances should be written as:
# conditions = [λ₀[2] =>-λ₀[1]]
conditions = [λ₀[2] =>-λ₀[1]]

σ₀ = transpose(aexp.mat)*λ₀

style = "RNF"

if style == "Graph"
    DP.res = DP.res.+1
else
    for in = 1:n_aut
        if style == "RNF"
            DP.res[in,:] = convert(Vector{Int},mysub(abs.(σ₀).-abs(λ₀[in]),conditions) .== 0)
        end
        if style == "CNF"
            DP.res[in,:] = convert(Vector{Int},mysub(σ₀.-λ₀[in],conditions) .== 0)
        end
    end
end


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                            Computation of parametrisation                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# pk denotes the order that is about to be computed
# nk the number of sets in that order
# ind_setk the iterative index from 1:nk
# ind_setG the global index relative to ind_setk


#~~~~~~~~~~~~~~~~~#
#           Order 0                  #
#~~~~~~~~~~~~~~~~~#
#
# the order is zero
p0 = 0;
# the number of sets in order zero is only one 
n0 = 1;
# the local to p0=0 index set (ind_set0) is then only one:
ind_set0 = 1 
# find equilibrium position with the function
# compute_order_zero
# that will find Y0 such that 
#    B.Y0 + (Q.Y0).Y0 + C0 = 0
# ATTENTION: there must be multiple solutions to this problem 
# the second entry of the function selects which solution to extract
# if C0 = 0, the function compute_order_zero returns automatically Y0 = 0
Y0 = compute_order_zero(sys,1)
# the global index_set relative to p0 = 0 and ind_set0 = 1) is found through
ind_setG0 = aexp.get(aexp.get([p0 ind_set0]))
# the computed Y0 is the mapping relating to the zero monomial
# so it is assigned to the DP.W at the entry ind_set:
DP.W[:,ind_setG0] = Y0
for i_full=1:n_full
    if DP.W[i_full,ind_setG0] != 0
        DP.Ws[i_full,ind_setG0] = Sym("W"*string(i_full)*"|"*string(ind_setG0))
        DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"|"*string(ind_setG0))=>DP.W[i_full,ind_setG0])]
    end
end
# if the function compute_order_zero struggles
# and if the right solution to the equation B.Y0 + (Q.Y0).Y0 + C0 = 0
# is known by the user, it can be directly input as:
#W0 = [W0₁ ...]

# the matrix B is now updated by adding to it 
# the gradient of the quadratic RHS calculated at Y0
# to do so, the function Q_fun is used
# Q_fun is a function that takes two vectors and a sparse tensor
# and computes Q_fun(Y1,Y2,Q) = (Q.Y1).Y2
# here the temporary single entry function ∇Q0_fun is defined
# which compute Q_fun at Y0 and generic Y
∇Q0_fun(Y) = Q_fun(Y,Y0,sys.Q)
# then the matrix ∇Q0 is extracted
∇Q0 = extract_Lin(∇Q0_fun,n_full)
# finally the linear RHS matrix B is updated by adding ∇Q0
sys.B+= ∇Q0
# ATTENTION: if the original B matrix is needed,
# it has to be saved before its update because it is no longer in sys


#~~~~~~~~~~~~~~~~~#
#           Order 1                  #
#~~~~~~~~~~~~~~~~~#
# AUTONOMOUS PART
#
# compute the full eigenvalues vector Λ, 
# right and left full eigenvectors matrices YR and YL
# Λ,YR,YL fulfil the conditions:
# (1)    sys.A*YR[:,k]*Λ[k] = sys.B*YR[:,k]
# (2)    Λ[k]*transpose(YL[:,k])*sys.A = transpose(YL[:,k])*sys.B
Λ,YR,YL = generalised_eigenproblem(sys)
# for the degrees of freedom of the full system 
# relating to the algebraic part of the DAE,
# the eigenvalues are set to the symbolic variable ∞
# the conditions (1) and (2) might not be respected for Λ[k]=∞
# before chosing which eigenvalues/vectors are the master,
# it is advised to simplify their expression
# in particular, in the non-conservative case, the simplification:
# sqrt(ξᵢ^2 - 1) => I*δᵢ            1 - ξᵢ^2 => δᵢ^2
# is recommended to keep the equations compact
# simpl = [Dict(sqrt(ξ[i]^2 - 1)=>im*δ[i]) for i=1:n_osc]
# Λ,YR,YL = mysub(Λ,simpl),mysub(YR,simpl),mysub(YL,simpl)
# 
# Normally the eigenvalues are organised as 
# the first n_aux eigenvalues Λ[1:n_aux] are equal to ∞
# then the meaningful eigenvalues starts usually sorted 
# by lower ω to higher
# Λ[n_aux+1:n_aux+2] are relative to ω₁
# Λ[n_aux+3:n_aux+4] are relative to ω₂
# and so on
# here the master are chosen as those relating to ω₁:
#yR = YR[:,n_aux+1:n_aux+2]
yR = YR[:,n_aux+2:-1:n_aux+1] #?? to test for having +i omega first....
#yL = YL[:,n_aux+1:n_aux+2]
yL = YL[:,n_aux+2:-1:n_aux+1] #?? to test
#test normalisation - a revoir // tester //OK this is the good one !
yR[:,1]=-yR[:,1]/(1/ω[1])/im
yR[:,2]=yR[:,2]/(1/ω[1])/im

#λ = Λ[n_aux+1:n_aux+2]
λ = Λ[n_aux+2:-1:n_aux+1] #good one here!
# any choice is possible but the sorting is not known before launching the script!
# one must then look at the eigenvalues sorting, then choose the masters
# after having chosen the master, 
# the linear part of the autonomous parametrisation can be filled
#
# the order is one
p1 = 1;
# ATTENTION: since the nonautonomous part is solved through
# an homological equation, here only the first n_aut indexes are filled
# the ind_set1 will then only iterate from 1 to n_aut
# leaving the remainder n_aut+1:n_rom for later
for ind_set1 = 1:n_aut
    # find the global ind_set relating to ind_set1
    ind_setG1 = aexp.get(aexp.get([p1 ind_set1]))
    # assign the chosen master right eigenvector to the corresponding DP.W
    DP.W[:,ind_setG1] = yR[:,ind_set1]
    for i_full=1:n_full
        if DP.W[i_full,ind_setG1] != 0
            DP.Ws[i_full,ind_setG1] = Sym("W"*string(i_full)*"|"*string(ind_setG1))
            DP.subs = [DP.subs;Dict(Sym("W"*string(i_full)*"|"*string(ind_setG1))=>DP.W[i_full,ind_setG1])]
        end
    end
    # assign the chosen master eigenvalue to the corresponding DP.f
    DP.f[ind_set1,ind_setG1] = λ[ind_set1]
    DP.fs[ind_set1,ind_setG1] = Sym("λ"*string(ind_set1))
    DP.subs = [DP.subs;Dict(Sym("λ"*string(ind_set1))=>λ[ind_set1])]
    # compute the matrix A*yR[aut]
    # which will be used for the top right border of the homological matrix
    yRs = 0*yR[:,ind_set1]
    for i_full=1:n_full
        if yR[i_full,ind_set1] != 0
            yRs[i_full] = Sym("yR"*string(i_full)*"|"*string(ind_set1))
            DP.subs = [DP.subs;Dict(Sym("yR"*string(i_full)*"|"*string(ind_set1))=>yR[i_full,ind_set1])]
        end
    end
    DP.AYR[:,ind_set1] = sys.A*yRs# yR[:,ind_set1]
    # compute the matrix yL[aut]ᵀ*A
    # which will be used for the bottom left border of the homological matrix
    yLsᵀ = 0*transpose(yL[:,ind_set1])
    for i_full=1:n_full
        if yL[i_full,ind_set1] != 0
            yLsᵀ[i_full] = Sym("yL"*string(i_full)*"|"*string(ind_set1))
            DP.subs = [DP.subs;Dict(Sym("yL"*string(i_full)*"|"*string(ind_set1))=>yL[i_full,ind_set1])]
        end
    end
    DP.YLᵀA[ind_set1,:] = yLsᵀ*sys.A    #transpose(yL[:,ind_set1])*sys.A    
end




#~~~~~~~~~~~~~~~~~#
#           Order 1                  #
#~~~~~~~~~~~~~~~~~#
# NON-AUTONOMOUS PART

if n_nonaut>0
    # augment λ with eigenvalues of the nonautonomous part:
    λ = [λ;im*symbols("Ω",positive=true);-im*symbols("Ω",positive=true)]
    λ = reshape(λ,1,n_rom)
    # assign the eigenvalues of the nonautonomous part to f:
    DP.f[n_aut+1,aexp.get(aexp.get([p1 n_aut+1]))] = λ[n_aut+1]
    DP.fs[n_aut+1,aexp.get(aexp.get([p1 n_aut+1]))] = Sym("λ"*string(n_aut+1))
    DP.subs = [DP.subs;Dict(Sym("λ"*string(n_aut+1))=>λ[n_aut+1])]
    DP.f[n_aut+2,aexp.get(aexp.get([p1 n_aut+2]))] = λ[n_aut+2]
    DP.fs[n_aut+2,aexp.get(aexp.get([p1 n_aut+2]))] = Sym("λ"*string(n_aut+2))
    DP.subs = [DP.subs;Dict(Sym("λ"*string(n_aut+2))=>λ[n_aut+2])]
    # assign the C⁺ₑₓₜ and C⁻ₑₓₜ to the RHS of the nonautonomous homological:
    DP.RHS_d[1:n_full,aexp.get(aexp.get([p1 n_aut+1]))] = C⁺ₑₓₜ
    DP.RHS_d[1:n_full,aexp.get(aexp.get([p1 n_aut+2]))] = C⁻ₑₓₜ
    # calculate σ and assign it to DP:
    DP.σ = λ*aexp.mat
    # σ is equal to:
    # [0;                                                   # order 0
    #  λ₁; λ₂; .. λₙᵣₒₘ;                                  # order 1
    #  2λ₁; λ₁+λ₂; .. 2λₙᵣₒₘ;                         # order 2
    #  3λ₁; 2λ₁+λ₂; .. 3λₙᵣₒₘ;                       # order 3
    #  ... ]                                                 # and so on
    # solve homological for the linear nonautonomous part:
    for ind_set1_nonaut = n_aut+1:n_aut+n_nonaut
        ind_setG1_nonaut = aexp.get(aexp.get([p1 ind_set1_nonaut]))
        println("solving order "*string(p1)*" and set "*string(ind_setG1_nonaut)*" with exponents:")
        println(aexp.mat[:,ind_setG1_nonaut])
        solve_homological!(ind_setG1_nonaut,DP,aexp,sys)
    end
else
    λ = reshape(λ,1,n_rom)
    # calculate σ and assign it to DP:
    DP.σ = λ*aexp.mat
end

#~~~~~~~~~~~~~~~~~#
#           Order p                  #
#~~~~~~~~~~~~~~~~~#
for p=2:o
    println("solving order "*string(p))
    fill_RHS_dyn!(p,DP,aexp,sys)
    fill_RHS_quad!(p,DP,aexp,sys)
    np=aexp.get(p)
    for ind_setp = 1:np
        ind_setGp = aexp.get(aexp.get([p ind_setp]))
        println("solving set "*string(ind_setGp)*" with exponents:")
        println(aexp.mat[:,ind_setGp])
        fill_RHS_lin!(aexp.mat[:,ind_setGp],DP,aexp,sys)
        solve_homological!(ind_setGp,DP,aexp,sys)
    end
end


#les lignes suivantes pour affichage de qq coefs individuels
#println("nonlinear term f such that ∂ₜz1 = [..] + f*z1^2*z2")
#println(mysub(mysub(mysub([DP.f[1,8]],DP.subs[end:-1:1]), [Dict(sqrt(ξ[i]^2 - 1)=>im*δ[i]) for i=1:n_osc]), [Dict(2*ξ[i]^3 - 2*ξ[i] =>-2*ξ[i]δ[i]^2) for i=1:n_osc]))

#println("nonlinear term f such that ∂ₜz2 = [..] + f*z1*z2^2")
#println(mysub(mysub(mysub([DP.f[2,9]],DP.subs[end:-1:1]), [Dict(sqrt(ξ[i]^2 - 1)=>im*δ[i]) for i=1:n_osc]), [Dict(2*ξ[i]^3 - 2*ξ[i] =>-2*ξ[i]δ[i]^2) for i=1:n_osc]))

#les lignes suivantes pour affichage des deux eqs complètes
z1=Sym("z_1");z2=Sym("z_2")
zₜ=sympy.zeros(n_aut,1)[:,1]
for i_ord=1:length(DP.f[1,:])
   for i_var=1:n_aut
       zₜ[i_var:i_var]+=mysub(mysub(mysub([DP.f[i_var,i_ord]],DP.subs[end:-1:1]), [Dict(sqrt(ξ[i]^2 - 1)=>im*δ[i]) for i=1:n_osc]), [Dict(2*ξ[i]^3 - 2*ξ[i] =>-2*ξ[i]δ[i]^2) for i=1:n_osc])*z1^aexp.mat[1,i_ord]*z2^aexp.mat[2,i_ord]
   end
end

#println("RHS dynamics whose LHS is zₜ=[∂z₁/∂t  ∂z₂/∂t  ... ]")
reduced_dynamics_latex_output(zₜ, "./test/Duffing_cubic_conservative_unforced_RNF_output.txt")
