using MAT
include("parametrisation.jl")
include("multiexponent.jl")

function matcont(DP::parametrisation_struct, aexp::multiexponent_struct)

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
            i_set = aexp.get(aexp.get([p i]))
            idof = 1
            # for idof=1:1 # For the moment only sdof systems are considered
                # dof=nodes[inode].dof[idof] 
                # if dof>0 
            mappings[index,1,idof] = DP.Wr[1,i_set]
            mappings_vel[index,1,idof] = DP.Wr[2,i_set]
                # end
            # end
            # for imode=1:info.neig # all computed modes
            #     mappings_modal_vel[index,imode]=V[:,imode]'*M*real.(P[p].Wr[1:info.neq,i])
            #     mappings_modal[index,imode]=V[:,imode]'*M*real.(P[p].Wr[info.neq+1:2*info.neq,i])
            # end
            index+=1
        end      
    end
    
    file = matopen("./output/param.mat","w")
    write(file, "ndof", 2)
    write(file, "nz", DP.n_aut)
    write(file, "max_order", DP.order)
    write(file, "mappings", mappings)
    # write(file, "mappings_modal",mappings_modal)
    write(file, "mappings_vel", mappings_vel)
    # write(file, "mappings_modal_vel",mappings_modal_vel)
    write(file, "Avector",Avector)
    write(file, "fdyn",fdyn)
    close(file)
    
end
    