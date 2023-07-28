using SymPy
#=using Combinatorics=#
using LinearAlgebra


function create_multiexp(n,o)
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
    mltexp=collect(multiexponents(n,0));
    MltExp=mltexp;
    a=Dict(mltexp[1]=> 1)
    a=merge(a, Dict([0 1]=>mltexp[1]))
    a=merge(a,Dict(0=>1))
    nc=1
    for p=1:o
        mltexp=collect(multiexponents(n,p));        
        MltExp=[MltExp;mltexp];
        a=merge(a,Dict(mltexp[i] => nc+i for i = 1:length(mltexp)))
        a=merge(a,Dict([p i] => mltexp[i] for i = 1:length(mltexp)))
        a=merge(a,Dict(p =>length(mltexp)))
        nc=nc+length(mltexp);
    end
    function get_a(a_set)
        return get(a,a_set,"-")
    end
    return get_a,a,MltExp,nc
end

function partitions_two(p)
    # partitions gives a list of all the 2 integers whose sum is p
    # expect the trivial one [p,0] and [0,p]
    #
    #
    # the part vector will be for example:
    #=
    
    partitions_two(5)        =
    4-element Vector{Vector{Int64}}:
    [4, 1]
    [3, 2]
    [2, 3]
    [1, 4]
    
    =#
    part=collect(partitions(p,2))
    if iseven(p)
        for k=length(part)-1:-1:1
            part=[part;[part[k][end:-1:1]]]
        end
    else
        for k=length(part):-1:1
            part=[part;[part[k][end:-1:1]]]
        end
    end
    return part
end

mutable struct multiexponent_struct
    get
    dic
    mat
    nc
end

function init_multiexponent_struct(n_rom,o)
    g,d,v,nc=create_multiexp(n_rom,o)
    return multiexponent_struct(g,d,v,nc)
end




function create_gen_vec(str,n)
    z=sympy.zeros(n,1)[:,1]
    for i=1:n
        z[i]=symbols(str*string(i))
    end
    return z
end






mutable struct system_struct
    A
    B
    Q
    C0
    # E
end

function extract_Quad(QuadFun,N::Int64)
    Q_tens=Dict{Any,Any};
    Y=create_gen_vec("Y",N)
    F=QuadFun(Y)
    for i=1:N
        for j=1:N
            for k=1:N
                g=diff(diff( F[i] ,Y[j]),Y[k]);
                if g !=0
                    Q_tens=merge(Q_tens,Dict([i,j,k]=>g/Sym(2)))
                end
            end
        end
    end
    return Q_tens
end

function extract_Lin(LinFun,N::Int64)
    Y=create_gen_vec("Y",N)
    F=LinFun(Y)
    L_mat=sympy.zeros(N)
    for i=1:N
        for j=1:N
            L_mat[i,j]=diff( F[i] ,Y[j]);
        end
    end
    return L_mat
end

function Q2_NL(Y1::Vector{Sym},Y2::Vector{Sym},sys::system_struct)
    F=sympy.zeros(n_full,1)[:,1];
    for i=1:n_full
        for j=1:n_full
            for k=1:n_full
                g=get(sys.Q,[i,j,k],Sym(0))
                F[i]+=g*Y1[j]*Y2[k]
            end
        end
    end
    return F
end





mutable struct parametrisation_struct
    RHS_d::Matrix{Sym}
    RHS_Q::Matrix{Sym}
    W::Matrix{Sym}    
    Ws::Matrix{Sym}
    fs::Matrix{Sym}
    f::Matrix{Sym}
end

function init_parametrisation_struct(n_full::Int64,n_rom::Int64,nc::Int64)
    mat_Nnc=sympy.zeros(n_full,nc);
    mat_nnc=sympy.zeros(n_rom,nc);
    return parametrisation_struct(mat_Nnc,mat_Nnc,mat_Nnc,mat_Nnc,mat_nnc,mat_nnc)
end




function fill_RHS_quad!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    part=partitions_two(p)
    # partitions gives a list of all the 2 integers whose sum is p
    # expect the trivial one [p,0] and [0,p]
    for k_part=1:length(part)
        p_L=part[k_part][1];
        p_R=part[k_part][2];
        n_L=aexp.get(p_L);
        n_R= aexp.get(p_R);
        for i_L=1:n_L
            for i_R=1:n_R
                set_L= aexp.get([p_L i_L])
                set_R=  aexp.get([p_R i_R])
                println(string(set_L)*" "*string(set_R))
                # the case of symmetric Q is not treated for generality
                #  RHS_quad(set_L+set_R) += Q(W_map(set_L), W_map(set_R))
                RHS_ind = aexp.get(set_L+set_R)
                psiL_ind = aexp.get(set_L)
                psiR_ind = aexp.get(set_R)
                DP.RHS_Q[:,RHS_ind] += Q2_NL(DP.W[:,psiL_ind],DP.W[:,psiR_ind],sys)
            end
        end
    end
    return nothing
end

function fill_RHS_dyn!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct)
    part=partitions_two(p+1)    
    # partitions gives a list of all the 2 integers whose sum is p
    # expect the trivial one [p+1,0] and [0,p+1]
    part=part[2:end-1][:]   
    # NB for forced this has to change! f^{s}_{s+nA} with  s \in nA is nonzero but still a nasty
    for k_part=1:length(part)        
        p_L=part[k_part][1];            println(p_L)
        p_R=part[k_part][2];
        n_L=aexp.get(p_L);
        n_R= aexp.get(p_R);
        for i_L=1:n_L
            for i_R=1:n_R
                set_L= aexp.get([p_L i_L])
                set_R=  aexp.get([p_R i_R])
                for s=1:length(set_L)
                    if set_L[s]>0
                        e_s=0*set_L;e_s[s]=1
                        # println(string(set_L)*" "*string(s)*" "*string(set_R)*" = "*string(set_L-e_s+set_R)*string(set_L[s]))
                        # 
                        # RHS_dyn(set_L-e_s+set_R) += - W(set_L)*f(s,set_R))*set_L[s]
                        RHS_ind = aexp.get(set_L-e_s+set_R)
                        W_ind = aexp.get(set_L)
                        f_ind = aexp.get(set_R)
                        DP.RHS_d[:,RHS_ind] += - DP.W[:,W_ind]*DP.f[s,f_ind]*set_L[s]
                    end
                end                
            end
        end
    end
    return nothing
end

function compute_order_zero(sys::system_struct,sol::Int64)
    N=size(sys.A)[1]
    if sys.C0==sympy.zeros(N,1)[:,1]
        x0=C0
    else
        x=create_gen_vec("x",N)
        x0_sol= solve(Q2_NL(x,x,sys)+sys.B*x+sys.C0)
        println("Static Solutions:")
        println(x0_sol)
        x0=sympy.zeros(N,1)[:,1];
        for in=1:N
            x0[in]=get(x0_sol[sol],x[in],"-")
        end
    end
    return x0
end


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#= system definition with quadratic recast and first order


M Uₜ  =    M U
M Vₜ  =  - K U - C V - G₁₁ U₁² - G₁₂ U₁₂ - G₂₂ U₂² - H₁₁₁ R₁ U₁ - H₁₁₂ R₁ U₂ - H₁₂₂ U₁ R₂ - H₂₂₂ U₂ R₂  - F0
    0    =    R₁ - U₁²
    0    =    R₂ - U₂²


Y[1] =  U₁
Y[2] =  U₂
Y[3] =  V₁
Y[4] =  V₂
Y[5] =  R₁
Y[6] =  R₂

The user can write the equations in the functions:
LHS_Lin(Yₜ) = RHS_Lin(Y) + RHS_Quad(Y) + C0

Upon extraction, the system will read:
A.Yₜ = B.Y + (Q.Y).Y + C0 

Where the matrices A and B, the tensor Q, and the vector C0
    will be stored in sys as: 
sys.A::Matrix{Sym}
sys.B::Matrix{Sym}
sys.Q::Dict{Any,Any}
sys.C0::Vector{Sym}


A  =[   [M [0 0 0 0;0 0 0 0]]
           [[0 0;0 0] M [0 0;0 0]]
           [0 0 0 0 0 0;0 0 0 0 0 0]  ]
B = [   [M [0 0 0 0;0 0 0 0]]
           [-C  -K [0 0;0 0]]
           [[0 0 0 0;0 0 0 0] I  ]
C0=[F0;0;0;0;0]

=#
M= Matrix{Sym}([[1 0];[0 1]])
ω₁=symbols("ω₁",positive=true)#integer=true, real=true
ω₂=symbols("ω₂",positive=true)#integer=true, real=true
K=[[ω₁^2 0];[0 ω₂^2]]
ωₙ=[ω₁;ω₂]
# create the diagonalised damping matrix either as
# generic diagonal damping:
ξ₁=symbols("ξ₁",positive=true)
ξ₂=symbols("ξ₂",positive=true)  
ξₙ=[ξ₁;ξ₂]
ζ=2*ξₙ.*ωₙ
C=diagm(ζ)
# definition of δₙ = √(1-ξₙ^2) used later for simplification
δₙ=sympy.zeros(size(M)[1])[:,1]
for i=1:size(M)[1]
        δₙ[i]=symbols("δ"*string(ξₙ[i]),positive=true)
end

n_rom=2
n_full=2*size(M)[1]+2
o=3

function LHS_Lin(Yₜ)
    F=sympy.zeros(n_full,1)[:,1];
    Uₜ=Yₜ[1:2]
    Vₜ=Yₜ[3:4]
    F[1:2]=M*Uₜ
    F[3:4]=M*Vₜ
    return F
end
function RHS_Lin(Y)
    F=sympy.zeros(n_full,1)[:,1];
    U=Y[1:2]
    V=Y[3:4]
    R=Y[5:6]
    F[1:2]=M*V
    F[3:4]=-C*V-K*U
    F[5:6]=R
    return F
end
function RHS_Quad(Y)
    F=sympy.zeros(n_full,1)[:,1];
    U₁=Y[1]
    U₂=Y[2]
    R₁=Y[5]
    R₂=Y[6]
    G¹₁₁=Sym("G¹₁₁")
    G²₁₁=Sym("G²₁₁");  G¹₁₂ = 2*G²₁₁;
    G¹₂₂=Sym("G¹₂₂");  G²₁₂ = 2*G¹₂₂;
    G²₂₂=Sym("G²₂₂")
    H¹₁₁₁=Sym("H¹₁₁₁")
    H²₁₁₁=Sym("H²₁₁₁"); H¹₁₁₂ = 3*H²₁₁₁;
    H¹₁₂₂=Sym("H¹₁₂₂"); H²₁₁₂ = H¹₁₂₂;
    H¹₂₂₂=Sym("H¹₂₂₂"); H²₁₂₂ = 3*H¹₂₂₂;
    H²₂₂₂=Sym("H²₂₂₂")
    F[1]=0;
    F[2]=0;
    F[3]=+        - (G¹₁₁*U₁^2   + G¹₁₂*U₂*U₁   + G¹₂₂*U₂^2) +
            -  (H¹₁₁₁*U₁*R₁ + H¹₁₁₂*U₂*R₁ + H¹₁₂₂*R₂*U₁ + H¹₂₂₂*R₂*U₂)
    F[4]=+        - (G²₁₁*U₁^2   + G²₁₂*U₂*U₁   + G²₂₂*U₂^2) +
            - (H²₁₁₁*U₁*R₁ + H²₁₁₂*U₂*R₁ + H²₁₂₂*R₂*U₁ + H²₂₂₂*R₂*U₂)
    F[5]=-U₁^2
    F[6]=-U₂^2
    return F
end
C0=sympy.zeros(n_full,1)[:,1]

sys=system_struct(extract_Lin(   LHS_Lin,    n_full),
                            extract_Lin(   RHS_Lin,    n_full),
                            extract_Quad(RHS_Quad,n_full),
                            C0)

aexp=init_multiexponent_struct(n_rom,o)

DP=init_parametrisation_struct(n_full,n_rom,aexp.nc)

p0=0;ic=1
DP.W[:,aexp.get(aexp.get([p0 ic]))]=compute_order_zero(sys,1)

#=
X=Sym("W^0_")
for k=1:2
multiexp=collect(multiexponents(n_rom,k))
for i=1:length(multiexp)
    X += symbols("W^"*string(k)*"_"*string(multiexp[i][1],multiexp[i][2],multiexp[i][3]))*prod(z.^multiexp[i])
end
end
=#

