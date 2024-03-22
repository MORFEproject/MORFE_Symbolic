function fill_RHS_quad!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    part = partitions_two(p)
    # partitions gives a list of all the 2 integers whose sum is p
    # expect the trivial one [p,0] and [0,p]
    # the trivial ones must be taken into account in the LHS
    # because they change B to B+∇Q_fun(W[0])
    for k_part = 1:length(part)
        p_L = part[k_part][1];
        p_R = part[k_part][2];
        n_L = aexp.get(p_L);
        n_R = aexp.get(p_R);
        for i_L = 1:n_L
            for i_R = 1:n_R
                set_L = aexp.get([p_L i_L])
                set_R = aexp.get([p_R i_R])
                # the case of symmetric Q is not treated for generality
                #  RHS_quad(set_L+set_R) += Q(W_map(set_L), W_map(set_R))
                ind_RHS = aexp.get(set_L+set_R)
                psiL_ind = aexp.get(set_L)
                psiR_ind = aexp.get(set_R)
                DP.RHS_Q[:,ind_RHS] += Q_fun(DP.Ws[:,psiL_ind],DP.Ws[:,psiR_ind],sys.Q)
            end
        end
    end
    return nothing
end

function fill_RHS_dyn!(p::Int64,DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    part = partitions_two(p+1)    
    # partitions gives a list of all the 2 integers whose sum is p+1
    # expect the trivial one [p+1,0] and [0,p+1]
    # those do not appear because f[0] = 0 and ∂ₜW[0] = 0
    # the outer ones [p,1] and [1,p]
    # are the ones that give:
    # 
    # LHS_dyn(set_current - e_s + e_r) = A*W(set_current)*f(s,r))*set_current[s]
    # in the autonomous case:
    # LHS_dyn(set_current) = A*W(set_current)*f(s,s))*set_current[s]
    # with σ = sum  f(s,s))*set_current[s]
    # 
    # and:
    # LHS_dyn(e_s - e_s + set_current) = A*W(e_s)*f(s,set_current))
    #
    # here I only take the inner ones:
    part = part[2:end-1][:]
    # NB for forced this has to change! 
    # f^{s}_{s+nA} with  s \in nA is nonzero but ...
    # the p-mapping W is still known from previous calculations of the automonous
    for k_part = 1:length(part)        
        p_L = part[k_part][1];# order of the mapping
        p_R = part[k_part][2];# order of the dynamics
        n_L = aexp.get(p_L);# number of sets in the mappings
        n_R = aexp.get(p_R);# number of sets in the dynamics
        for i_L = 1:n_L
            for i_R = 1:n_R
                set_L = aexp.get([p_L i_L])# get the set for the mappings
                set_R = aexp.get([p_R i_R])# get the set for the dynamics
                for s = 1:length(set_L)       # for s = 1:n_rom
                    if set_L[s]>0               # if the exponent of zₛ is nonzero
                        e_s = 0*set_L;e_s[s] = 1    # create unit vector e_s
                        # the term of the mapping:
                        # W(set_L)*z^(set_L)
                        # when derived in time (along zₛ) gives a term:
                        # W(set_L)*f(s,set_R))*set_L[s]*z^(set_L-e_s+set_R)
                        # then it will go on the RHS with a minus and premult by A:
                        # RHS_dyn(set_L-e_s+set_R) += - A*W(set_L)*f(s,set_R))*set_L[s]
                        ind_f = aexp.get(set_R)
                        fs_set_R = DP.fs[s,ind_f]
                        #if fs_set_R != 0    # only fill if the reduced dyn is nonzero
                            ind_RHS = aexp.get(set_L-e_s+set_R)
                            ind_W = aexp.get(set_L)
                            DP.RHS_d[:,ind_RHS] += - sys.B*DP.Ws[:,ind_W]*fs_set_R*set_L[s]
                        #end
                    end
                end                
            end
        end
    end
    return nothing
end

function fill_RHS_lin!(set_current::Vector{Int},DP::parametrisation_struct,aexp::multiexponent_struct,sys::system_struct)
    # 
    # LHS_dyn(set_dependent - e_s + e_r) = A*W(set_dependent)*f(s,r))*set_dependent[s]
    # sum over s ∈ set_dependent, sum over r ∈ (n_aut+1,n_rom)
    # f(s,r) = δₛᵣ λᵣ    if s and r ∈ (1,n_aut)
    # f(s,r)  != 0      if s ∈ (1,n_aut) and r ∈ (n_aut+1,n_rom)
    # in the autonomous case:
    # LHS_dyn(set_current) = A*W(set_current)*f(s,s))*set_current[s]
    # with σ = sum  f(s,s))*set_current[s]
    # 
    # in the nonautonomous case:
    # for r ∈ n_nonaut
    #   if r ∈ set_current
    #       for s ∈ n_aut
    #           set_dependent = set_current + e_s - e_r   
    #           LHS_dyn(set_current) = A*W(set_dependent)*f(s,r))*set_dependent[s]
    #       end
    #   end
    # end
    ind_RHS = aexp.get(set_current)
    for r = 1+DP.n_aut:DP.n_rom
        if set_current[r]>0
            e_r = 0*set_current;e_r[r] = 1    # create unit vector e_r
            ind_f = aexp.get(e_r)
            for s = 1:DP.n_aut                
                fs_r = DP.f[s,ind_f]                    
                if fs_r != 0    # only fill if the reduced dyn is nonzero
                    e_s = 0*set_current;e_s[s] = 1    # create unit vector e_s
                    set_dependent = set_current + e_s - e_r
                    ind_W = aexp.get(set_dependent)
                    DP.RHS_d[:,ind_RHS] += - sys.B*DP.Ws[:,ind_W]*fs_r*set_dependent[s]
                end
            end                
        end
    end
    return nothing
end