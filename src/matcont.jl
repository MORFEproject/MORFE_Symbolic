function matcont(info::Sinfo,nodes::Vector{Snode},M::SparseMatrixCSC{Float64},V::Matrix{Float64},P::Vector{Parametrisation})

    howmany = 0
    for p in 0:DP.order   
        howmany += aexp.get(p)
    end
    
    Avector = zeros(howmany,DP.n_rom)
    fdyn = zeros(howmany,DP.n_aut)
    mappings = zeros(howmany,1,3) # At present restricted to systems with only one "node"
    mappings_vel = zeros(howmany,1,3)
    # mappings_modal = zeros(howmany,info.neig)
    # mappings_modal_vel = zeros(howmany,info.neig)
    
    index = 1
    for p in 0:DP.order   
        for i in 1:aexp.get(p)
            Avector[index,:] = aexp.get([p i])
            index += 1
        end
    end
    
    index = 1
    for p in 0:DP.order  
        for i in 1:aexp.get(p)
            i_set = aexp.get(aexp.get([p i]))
            for j in 1:DP.n_aut
                fdyn[index,j] = DP.fr[j,i_set]
            end
            index += 1
        end
    end
    
    index = 1
    for p in 0:DP.order 
        for i in 1:aexp.get(p)
            for idof=1:2 
                dof=nodes[inode].dof[idof] 
                if dof>0 
                    mappings_vel[index,inode,idof]=real(P[p].Wr[dof,i])
                    mappings[index,inode,idof]=real(P[p].Wr[info.neq+dof,i])
                end
            end
            for imode=1:info.neig # all computed modes
                mappings_modal_vel[index,imode]=V[:,imode]'*M*real.(P[p].Wr[1:info.neq,i])
                mappings_modal[index,imode]=V[:,imode]'*M*real.(P[p].Wr[info.neq+1:2*info.neq,i])
            end
            index+=1
        end      
    end
    
    file = matopen("./output/param.mat","w")
    write(file, "ndof",2*info.neq)
    write(file, "nz",info.nza)
    write(file, "max_order",info.max_order)
    write(file, "mappings",mappings)
    write(file, "mappings_modal",mappings_modal)
    write(file, "mappings_vel",mappings)
    write(file, "mappings_modal_vel",mappings_modal_vel)
    write(file, "Avector",Avector)
    write(file, "fdyn",fdyn)
    close(file)
    
end
    