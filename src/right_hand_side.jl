"""
In order to understand the functions in this module, references to the following paper will be made:
A. Vizzaccaro, G. Gobat, A. Frangi, C. Touze'. Direct parametrisation of invariant manifolds 
for non-autonomous forced systems including superharmonic resonances (2023).
"""

"""
Function to fill the part of the RHS of the homological equations related to [Q(W,W)]ₚ.
Takes as input arguments:\\
    - p: The current order of parametrisation that is being solved for;\\
    - DP: The parametrisation structure;\\
    - aexp: The multiexponent structure;\\
    - sys: The system of equations structure.\\
Does not return anything, but alters DP.
"""
function fill_RHS_quad!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    part = partitions_two(p)
    for k_part = 1:length(part)
        p_L = part[k_part][1];                 # order of the left mapping
        p_R = part[k_part][2];                 # order of the right mapping 
        n_L = aexp.get(p_L);                   # number of sets in the left mapping
        n_R = aexp.get(p_R);                   # number of sets in the right mapping
        for i_L = 1:n_L
            for i_R = 1:n_R
                set_L = aexp.get([p_L i_L])
                set_R = aexp.get([p_R i_R])
                # the case of symmetric Q is not treated for generality
                ind_RHS = aexp.get(set_L+set_R)
                psiL_ind = aexp.get(set_L)
                psiR_ind = aexp.get(set_R)
                DP.RHS_Q[:,ind_RHS] += Q_fun(DP.Ws[:,psiL_ind],DP.Ws[:,psiR_ind],sys.Q)
            end
        end
    end
    return nothing
end

"""
Function to fill the part of the RHS of the homological equations related to N₃.
Takes as input arguments:\\
    - p: The current order of parametrisation that is being solved for;\\
    - DP: The parametrisation structure;\\
    - aexp: The multiexponent structure;\\
    - sys: The system of equations structure.\\
Does not return anything, but alters DP.
"""
function fill_RHS_dyn!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    part = partitions_two(p+1)    
    part = part[2:end-1][:]
    for k_part = 1:length(part)        
        p_L = part[k_part][1];                      # order of the mapping
        p_R = part[k_part][2];                      # order of the dynamics
        n_L = aexp.get(p_L);                        # number of sets in the mapping
        n_R = aexp.get(p_R);                        # number of sets in the dynamics
        for i_L = 1:n_L
            for i_R = 1:n_R
                set_L = aexp.get([p_L i_L])
                set_R = aexp.get([p_R i_R])
                for s = 1:length(set_L)       
                    if set_L[s]>0                   # if the exponent of zₛ is nonzero
                        e_s = 0*set_L;e_s[s] = 1    # create unit vector e_s
                        ind_f = aexp.get(set_R)
                        fs_set_R = DP.fs[s,ind_f]
                        if fs_set_R != 0            # only fill if the reduced dyn is nonzero
                            ind_RHS = aexp.get(set_L-e_s+set_R)
                            ind_W = aexp.get(set_L)
                            DP.RHS_d[:,ind_RHS] += - sys.B*DP.Ws[:,ind_W]*fs_set_R*set_L[s]
                        end
                    end
                end                
            end
        end
    end
    return nothing
end

"""
Function to fill the part of the RHS of the homological equations related to N₂.
Takes as input arguments:\\
    - set_current: Corresponds to α(p,k) in the notation of the paper cited at
    the beginning of the module;\\
    - DP: The parametrisation structure;\\
    - aexp: The multiexponent structure;\\
    - sys: The system of equations structure.\\
Does not return anything, but alters DP.
"""
function fill_RHS_lin!(set_current::Vector{Int},DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    ind_RHS = aexp.get(set_current)
    for r = 1+DP.n_aut:DP.n_rom
        if set_current[r]>0                                     # if the exponent of zᵣ is nonzero
            e_r = 0*set_current;e_r[r] = 1                      # create unit vector e_r
            ind_f = aexp.get(e_r)
            for s = 1:DP.n_aut                
                fs_r = DP.f[s,ind_f]                    
                if fs_r != 0                                    # only fill if the reduced dyn is nonzero
                    e_s = 0*set_current;e_s[s] = 1              # create unit vector e_s
                    set_dependent = set_current + e_s - e_r     # α(pW,k) in the notation of the reference paper
                    ind_W = aexp.get(set_dependent)
                    DP.RHS_d[:,ind_RHS] += - sys.B*DP.Ws[:,ind_W]*fs_r*set_dependent[s]
                end
            end                
        end
    end
    return nothing
end