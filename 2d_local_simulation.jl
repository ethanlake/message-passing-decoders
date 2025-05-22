include("common_functions.jl")

"""
simulation of the 2D local message passing decoder
"""

function update_fields!(fields,new_fields,anyons,bconds,onenorm)
    """
    fields: LxLx2x2 of Ints (fields in each direction). last indices are (x/y, positive/negative)
    """
    L = size(fields)[1]
    ind(i) = mod1(i,L)

    rmin = bconds == "periodic" ? 1 : 2
    rmax = bconds == "periodic" ? L : L-1

    for i in rmin:rmax, j in rmin:rmax
        ip1 = ind(i+1); im1 = ind(i-1)
        jp1 = ind(j+1); jm1 = ind(j-1)

        if onenorm 

            #### 1-norm distance fields #### 
            # currently the lightfronts propagate with the 1-norm, but they are erased with the \infty norm; getting them to be erased with the 1-norm necessesitates a few extra variables...  

            ### +x fields ### 
            newfield = Inf 
            if anyons[im1,j] 
                newfield = 1 
            end 
            if fields[im1,j,1,1] != 0  
                newfield = min(newfield,fields[im1,j,1,1] + 1)
            end 
            if (fields[im1,jm1,1,1] != 0 || anyons[im1,jm1]) && fields[i,jm1,1,1] != 0 
                if anyons[im1,jm1] 
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[im1,jm1,1,1]+2)#,fields[i,jm1,1,1]+1)
                end 
            end 
            if (fields[im1,jp1,1,1] != 0 || anyons[im1,jp1]) && fields[i,jp1,1,1] != 0 
                if anyons[im1,jp1]
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[im1,jp1,1,1]+2)#,fields[i,jp1,1,1]+1) 
                end
            end 
            new_fields[i,j,1,1] = newfield == Inf ? 0 : newfield

            ### -x fields ###
            newfield = Inf 
            if anyons[ip1,j] 
                newfield = 1 
            end 
            if fields[ip1,j,1,2] != 0  
                newfield = min(newfield,fields[ip1,j,1,2] + 1)
            end 
            if (fields[ip1,jm1,1,2] != 0 || anyons[ip1,jm1]) && fields[i,jm1,1,2] != 0 
                if anyons[ip1,jm1] 
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[ip1,jm1,1,2]+2)#,fields[i,jm1,1,2]+1)
                end 
            end 
            if (fields[ip1,jp1,1,2] != 0 || anyons[ip1,jp1]) && fields[i,jp1,1,2] != 0 
                if anyons[ip1,jp1]
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[ip1,jp1,1,2]+2)#,fields[i,jp1,1,2]+1) 
                end
            end 
            new_fields[i,j,1,2] = newfield == Inf ? 0 : newfield

            ### +y fields ###
            newfield = Inf
            if anyons[i,jm1] 
                newfield = 1 
            end
            if fields[i,jm1,2,1] != 0 
                newfield = min(newfield,fields[i,jm1,2,1] + 1)
            end
            if (fields[im1,jm1,2,1] != 0 || anyons[im1,jm1]) && fields[im1,j,2,1] != 0 
                if anyons[im1,jm1]
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[im1,jm1,2,1]+2)#,fields[im1,j,2,1]+1)
                end 
            end
            if (fields[ip1,jm1,2,1] != 0 || anyons[ip1,jm1]) && fields[ip1,j,2,1] != 0 
                if anyons[ip1,jm1]
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[ip1,jm1,2,1]+2)#,fields[ip1,j,2,1]+1)
                end 
            end
            new_fields[i,j,2,1] = newfield == Inf ? 0 : newfield

            ### -y fields ###
            newfield = Inf
            if anyons[i,jp1] 
                newfield = 1 
            end
            if fields[i,jp1,2,2] != 0 
                newfield = min(newfield,fields[i,jp1,2,2] + 1)
            end
            if (fields[im1,jp1,2,2] != 0 || anyons[im1,jp1]) && fields[im1,j,2,2] != 0 
                if anyons[im1,jp1]
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[im1,jp1,2,2]+2)#,fields[im1,j,2,2]+1)
                end 
            end
            if (fields[ip1,jp1,2,2] != 0 || anyons[ip1,jp1]) && fields[ip1,j,2,2] != 0 
                if anyons[ip1,jp1]
                    newfield = min(newfield,2)
                else 
                    newfield = min(newfield,fields[ip1,jp1,2,2]+2)#,fields[ip1,j,2,2]+1)
                end 
            end
            new_fields[i,j,2,2] = newfield == Inf ? 0 : newfield
        
        else 
        ### infinity norm ### 
            ### +x fields ####
            if anyons[im1,j] || anyons[im1,jm1] || anyons[im1,jp1]
                new_fields[i,j,1,1] = 1
            else 
                vals = [fields[im1,j,1,1],fields[im1,jm1,1,1],fields[im1,jp1,1,1]]
                nonzerovals = vals[vals .> 0]
                if length(nonzerovals) > 0
                    new_fields[i,j,1,1] = minimum(nonzerovals) + 1
                else 
                    new_fields[i,j,1,1] = 0
                end
            end

            ### -x fields ####
            if anyons[ip1,j] || anyons[ip1,jm1] || anyons[ip1,jp1]
                new_fields[i,j,1,2] = 1
            else 
                vals = [fields[ip1,j,1,2],fields[ip1,jm1,1,2],fields[ip1,jp1,1,2]]
                nonzerovals = vals[vals .> 0]
                if length(nonzerovals) > 0
                    new_fields[i,j,1,2] = minimum(nonzerovals) + 1
                else 
                    new_fields[i,j,1,2] = 0
                end
            end

            ### +y fields ####
            if anyons[i,jm1] || anyons[im1,jm1] || anyons[ip1,jm1]
                new_fields[i,j,2,1] = 1
            else 
                vals = [fields[i,jm1,2,1],fields[im1,jm1,2,1],fields[ip1,jm1,2,1]]
                nonzerovals = vals[vals .> 0]
                if length(nonzerovals) > 0
                    new_fields[i,j,2,1] = minimum(nonzerovals) + 1
                else 
                    new_fields[i,j,2,1] = 0
                end
            end

            ### -y fields ####
            if anyons[i,jp1] || anyons[im1,jp1] || anyons[ip1,jp1]
                new_fields[i,j,2,2] = 1
            else 
                vals = [fields[i,jp1,2,2],fields[im1,jp1,2,2],fields[ip1,jp1,2,2]]
                nonzerovals = vals[vals .> 0]
                if length(nonzerovals) > 0
                    new_fields[i,j,2,2] = minimum(nonzerovals) + 1
                else 
                    new_fields[i,j,2,2] = 0
                end
            end
            
        end 
    end 
    fields .= new_fields
end 

function anyons_source_fields!(anyons,fields) # ensures that fields are always updated in the 1-balls around the each anyon's position
    L = size(fields)[1]
    ind(i) = mod1(i,L)
    for i in 1:L, j in 1:L 
        if anyons[i,j]
            ip1 = ind(i+1); im1 = ind(i-1)
            jp1 = ind(j+1); jm1 = ind(j-1)
            fields[ip1,j,1,1] = 1 
            fields[ip1,jm1,1,1] = 1 
            fields[ip1,jp1,1,1] = 1 
            fields[im1,j,1,2] = 1
            fields[im1,jm1,1,2] = 1
            fields[im1,jp1,1,2] = 1
            fields[i,jp1,2,1] = 1
            fields[im1,jp1,2,1] = 1
            fields[ip1,jp1,2,1] = 1
            fields[i,jm1,2,2] = 1
            fields[im1,jm1,2,2] = 1
            fields[ip1,jm1,2,2] = 1
        end 
    end 

end 

function update_state!(state,correction,anyons,fields,new_fields,r,p,synch,bconds)
    """
    state, correction: LxLx2 of Bools 
    anyons: LxL of Bools
    fields, new_fields: LxLx2x2 of Ints (fields in each direction)
    r: ratio of field updates to spin updates
    p: error probability 
    synch: if true, updates are done synchronously 
    bconds: boundary conditions ∈ {"periodic","open"}. if open, 
    """

    L = size(fields)[1]
    ind(i) = mod1(i,L)

    rmin = bconds == "periodic" ? 1 : 2
    rmax = bconds == "periodic" ? L : L-1

    onenorm = ~true 

    if synch 
        for _ in 1:(r-1)
            update_fields!(fields,new_fields,anyons,bconds,onenorm)
        end 

        # update the anyon positions 
        correction .= false 
        for i in rmin:rmax, j in rmin:rmax  # 
            if anyons[i,j] 
                upd = bconds == "open" ? true : rand() < .9 # small amount of stochasticity added to break out of limit cycles at small system sizes with periodic boundary conditions
                if upd 
                    if maximum(fields[i,j,:,:]) != 0 # move somewhere 
                        im1 = ind(i-1); jm1 = ind(j-1)
                        mindist = minimum(fields[i,j,:,:][fields[i,j,:,:] .> 0])
                        if fields[i,j,1,1] == mindist 
                            correction[im1,j,1] = true 
                        elseif fields[i,j,2,1] == mindist 
                            correction[i,jm1,2] = true
                        elseif fields[i,j,2,2] == mindist 
                            correction[i,j,2] = true
                        elseif fields[i,j,1,2] == mindist 
                            correction[i,j,1] = true
                        end
                    end 
                end 
            end
        end 
        state .⊻= correction

        # add noise 
        for i in rmin:rmax, j in rmin:rmax, o in 1:2
            if rand() < p 
                state[i,j,o] ⊻= true 
            end 
        end

        # compute anyon positions (not efficient but fine for now)
        anyons .= get_synds(state)

        anyons_source_fields!(anyons,fields)

        update_fields!(fields,new_fields,anyons,bconds,onenorm)


    ### asychronous updates ### 
    else 
        function locally_source_fields!(x,y)
            fields[ind(x+1),y,1,1] = 1
            fields[ind(x-1),y,1,2] = 1 
            fields[x,ind(y+1),2,1] = 1
            fields[x,ind(y-1),2,2] = 1
        end 
        for _ in 1:L^2*r 
            i = rand(rmin:rmax); j = rand(rmin:rmax)
            ip1 = ind(i+1); im1 = ind(i-1)
            jp1 = ind(j+1); jm1 = ind(j-1)

            if rand() > 1/r # update fields 

                ### infinity norm ### 
                if ~onenorm 
                    ### +x fields ####
                    if anyons[im1,j] || anyons[im1,jm1] || anyons[im1,jp1]
                        fields[i,j,1,1] = 1
                    else 
                        vals = [fields[im1,j,1,1],fields[im1,jm1,1,1],fields[im1,jp1,1,1]]
                        nonzerovals = vals[vals .> 0]
                        if length(nonzerovals) > 0
                            fields[i,j,1,1] = minimum(nonzerovals) + 1
                        else 
                            fields[i,j,1,1] = 0
                        end
                    end

                    ### -x fields ####
                    if anyons[ip1,j] || anyons[ip1,jm1] || anyons[ip1,jp1]
                        fields[i,j,1,2] = 1
                    else 
                        vals = [fields[ip1,j,1,2],fields[ip1,jm1,1,2],fields[ip1,jp1,1,2]]
                        nonzerovals = vals[vals .> 0]
                        if length(nonzerovals) > 0
                            fields[i,j,1,2] = minimum(nonzerovals) + 1
                        else 
                            fields[i,j,1,2] = 0
                        end
                    end

                    ### +y fields ####
                    if anyons[i,jm1] || anyons[im1,jm1] || anyons[ip1,jm1]
                        fields[i,j,2,1] = 1
                    else 
                        vals = [fields[i,jm1,2,1],fields[im1,jm1,2,1],fields[ip1,jm1,2,1]]
                        nonzerovals = vals[vals .> 0]
                        if length(nonzerovals) > 0
                            fields[i,j,2,1] = minimum(nonzerovals) + 1
                        else 
                            fields[i,j,2,1] = 0
                        end
                    end

                    ### -y fields ####
                    if anyons[i,jp1] || anyons[im1,jp1] || anyons[ip1,jp1]
                        fields[i,j,2,2] = 1
                    else 
                        vals = [fields[i,jp1,2,2],fields[im1,jp1,2,2],fields[ip1,jp1,2,2]]
                        nonzerovals = vals[vals .> 0]
                        if length(nonzerovals) > 0
                            fields[i,j,2,2] = minimum(nonzerovals) + 1
                        else 
                            fields[i,j,2,2] = 0
                        end
                    end

                ### 1 norm ### 
                else 
                    ### +x fields ### 
                    newfield = Inf 
                    if anyons[im1,j] 
                        newfield = 1 
                    end 
                    if fields[im1,j,1,1] != 0  
                        newfield = min(newfield,fields[im1,j,1,1] + 1)
                    end 
                    if (fields[im1,jm1,1,1] != 0 || anyons[im1,jm1]) && fields[i,jm1,1,1] != 0 
                        if anyons[im1,jm1] 
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[im1,jm1,1,1]+2)#,fields[i,jm1,1,1]+1)
                        end 
                    end 
                    if (fields[im1,jp1,1,1] != 0 || anyons[im1,jp1]) && fields[i,jp1,1,1] != 0 
                        if anyons[im1,jp1]
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[im1,jp1,1,1]+2)#,fields[i,jp1,1,1]+1) 
                        end
                    end 
                    fields[i,j,1,1] = newfield == Inf ? 0 : newfield

                    ### -x fields ###
                    newfield = Inf 
                    if anyons[ip1,j] 
                        newfield = 1 
                    end 
                    if fields[ip1,j,1,2] != 0  
                        newfield = min(newfield,fields[ip1,j,1,2] + 1)
                    end 
                    if (fields[ip1,jm1,1,2] != 0 || anyons[ip1,jm1]) && fields[i,jm1,1,2] != 0 
                        if anyons[ip1,jm1] 
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[ip1,jm1,1,2]+2)#,fields[i,jm1,1,2]+1)
                        end 
                    end 
                    if (fields[ip1,jp1,1,2] != 0 || anyons[ip1,jp1]) && fields[i,jp1,1,2] != 0 
                        if anyons[ip1,jp1]
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[ip1,jp1,1,2]+2)#,fields[i,jp1,1,2]+1) 
                        end
                    end 
                    fields[i,j,1,2] = newfield == Inf ? 0 : newfield

                    ### +y fields ###
                    newfield = Inf
                    if anyons[i,jm1] 
                        newfield = 1 
                    end
                    if fields[i,jm1,2,1] != 0 
                        newfield = min(newfield,fields[i,jm1,2,1] + 1)
                    end
                    if (fields[im1,jm1,2,1] != 0 || anyons[im1,jm1]) && fields[im1,j,2,1] != 0 
                        if anyons[im1,jm1]
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[im1,jm1,2,1]+2)#,fields[im1,j,2,1]+1)
                        end 
                    end
                    if (fields[ip1,jm1,2,1] != 0 || anyons[ip1,jm1]) && fields[ip1,j,2,1] != 0 
                        if anyons[ip1,jm1]
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[ip1,jm1,2,1]+2)#,fields[ip1,j,2,1]+1)
                        end 
                    end
                    fields[i,j,2,1] = newfield == Inf ? 0 : newfield

                    ### -y fields ###
                    newfield = Inf
                    if anyons[i,jp1] 
                        newfield = 1 
                    end
                    if fields[i,jp1,2,2] != 0 
                        newfield = min(newfield,fields[i,jp1,2,2] + 1)
                    end
                    if (fields[im1,jp1,2,2] != 0 || anyons[im1,jp1]) && fields[im1,j,2,2] != 0 
                        if anyons[im1,jp1]
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[im1,jp1,2,2]+2)#,fields[im1,j,2,2]+1)
                        end 
                    end
                    if (fields[ip1,jp1,2,2] != 0 || anyons[ip1,jp1]) && fields[ip1,j,2,2] != 0 
                        if anyons[ip1,jp1]
                            newfield = min(newfield,2)
                        else 
                            newfield = min(newfield,fields[ip1,jp1,2,2]+2)#,fields[ip1,j,2,2]+1)
                        end 
                    end
                    fields[i,j,2,2] = newfield == Inf ? 0 : newfield
                end 


            else # update anyons 
                if anyons[i,j] 
                    if maximum(fields[i,j,:,:]) != 0 # move somewhere 
                        im1 = ind(i-1); jm1 = ind(j-1)
                        mindist = minimum(fields[i,j,:,:][fields[i,j,:,:] .> 0])
                        if fields[i,j,1,1] == mindist
                            state[im1,j,1] ⊻= true 
                            anyons[i,j] ⊻= true; anyons[im1,j] ⊻= true
                            if anyons[im1,j]
                                locally_source_fields!(im1,j)
                            end
                        elseif fields[i,j,2,1] == mindist
                            state[i,jm1,2] ⊻= true
                            anyons[i,j] ⊻= true; anyons[i,jm1] ⊻= true
                            if anyons[i,jm1]
                                locally_source_fields!(i,jm1)
                            end
                        elseif fields[i,j,2,2] == mindist 
                            state[i,j,2] ⊻= true
                            anyons[i,j] ⊻= true; anyons[i,ind(j+1)] ⊻= true
                            if anyons[i,ind(j+1)]
                                locally_source_fields!(i,ind(j+1))
                            end
                        elseif fields[i,j,1,2] == mindist 
                            state[i,j,1] ⊻= true
                            anyons[i,j] ⊻= true; anyons[ind(i+1),j] ⊻= true
                            if anyons[ind(i+1),j]
                                locally_source_fields!(ind(i+1),j)
                            end
                        end
                    end 
                end

                # add noise (noise frequency matched to anyon updates)
                if rand() < p 
                    state[i,j,1] ⊻= true 
                    anyons[i,j] ⊻= true; anyons[ind(i+1),j] ⊻= true
                end 
                if rand() < p 
                    state[i,j,2] ⊻= true 
                    anyons[i,j] ⊻= true; anyons[i,ind(j+1)] ⊻= true
                end 
            end 
        end
    end 
end 

function parameter_repository(mode,L,p,synch,vary_L)
    pmin = 0; pmax = 1; nps = 1; nLs = 1 
    samps = 1; samps_vec = [samps]; ps = [p]; Ls = [L]; Ts = [1]

    if mode == "trel"
        ### fix system size, vary p ###
        if ~vary_L 
            nps = 5
            samps = round(Int,5000 / L)
            samps_vec = [samps for _ in 1:nps]
            
            pmin = .0025; pmax = .04 
            ps = [10^(x) for x in LinRange(log10(pmin),log10(pmax),nps)] 

        ### vary system size, fix p ###
        else 
            Lmin = 8; Lmax = 128
            nels = 10
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nels)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,800000/el) for el in Ls] # poosamps 
        end
        Ts = [1 for _ in 1:nps]
    end 

    if mode == "Ft"
        nps = 11
        samps_vec = [round(Int,200000 / L) for _ in 1:nps]

        if ~vary_L  # vary p 
            Ts = [L for _ in 1:nps]
            pmin = .01; pmax = .1 
            ps = [p for p in LinRange(pmin,pmax,nps)]

        else # vary L  
            Lmin = 32; Lmax = 512
            nels = 10
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nels)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,800000/el) for el in Ls] 
            Ts = [el for _ in Ls]
        end 
    end 

    if mode == "erode"
        
        if ~vary_L 
            nps = 12
            samps_vec = [round(Int,355000 / L) for _ in 1:nps]
            if synch 
                pmin = .04; pmax = .1 
            else 
                pmin = .03; pmax = .07  # for infinity norm 
            end 
            ps = [p for p in LinRange(pmin,pmax,nps)]

        else # varying L 
            # for looking at the erosion time: can go to very large L 
            Lmin = 32; Lmax = 512
            nps = 7
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nps)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,305200/el) for el in Ls] # poosamps 

            # for getting the logical error rate: use smaller system sizes 
            Lmin = 8; Lmax = 64
            if p == .05 
                Lmax = 50 
            elseif p == .045 
                Lmax = Lmax = 35 
            end 
            if ~synch 
                if p == .04 
                    Lmax = 65
                elseif p == .035
                    Lmax = 55
                elseif p == .03 
                    Lmax = 45 
                end 
            end 
            nps = 5
            Ls = [round(Int,el) for el in LinRange(Lmin,Lmax,nps)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,7000000/el) for el in Ls] # poosamps # professional 
            samps_vec = [round(Int,700000/el) for el in Ls] # poosamps 

        end 
        Ts = [1 for _ in 1:nps]
    end        

    if mode ∈ ["erode_times" "hist"]
        ps = [pmin]; nps = 1  
    else
        if vary_L 
            ps = [p for _ in 1:length(Ls)]
        else 
            Ls = [L for _ in 1:nps]
        end 
    end 

    params = Dict{String, Any}()
    params["samps"] = samps_vec; params["Ts"] = Ts; params["mode"] = mode; params["L"] = L; params["nps"] = nps; params["ps"] = ps; params["Ls"] = Ls; params["p"] = p; params["vary_L"] = vary_L

    return params 
end 

function main()

    """
    supported modes: 
    * "trel": compute relaxation time for online decoding 
    * "Ft": get decoding fidelity and time-dependent magnetization for online decoding 
    * "hist": get history of particular evolution
    * "erode": compute erosion fidelities and erosion times (and statistics thereof) 
    """

    mode = "hist" # "trel" "Ft" "hist" "erode" 
    L = 75 # 
    p = .24 # (used if varying L) 
    vary_L = ~true # if true, vary system size; if false, use fixed system size and vary p  

    r = 3 # number of field updates per spin update. need r > 2 to rigorously guarantee linear erosion
    synch = true 
    bconds = "periodic" # "periodic" "open" 

    out_adj = ""

    params = parameter_repository(mode,L,p,synch,vary_L)

    Ts = params["Ts"]; samps = params["samps"]; ps = params["ps"]; nps = params["nps"]; p = params["p"];  Ls = params["Ls"]
    params["synch"] = synch; params["r"] = r; params["bconds"] = bconds

    data_keys = ["Ft" "Mt"] ∪ ["hist" "field_hist" "anyon_hist"] ∪ ["trels" "trel_stats"] ∪ ["erode_times" "erode_stats" "anyon_density"] 
    data = Dict{String, Any}(key => 0 for key in data_keys)

    println("details of simulation: ")
    println("synch = $synch")
    if vary_L 
        println("p = $p")
        println("Ls = $Ls")
    else 
        println("system size: $L")
        println("ps = $ps")
    end 
    println("mode = $mode")
    println("field update speed = $r")


    ### write history of evolution ### 
    if mode == "hist" 
        T = round(Int,L) 
    
        println("running for time T = $T")
        data["hist"] = zeros(Bool,T,L,L,2)
        data["field_hist"] = zeros(Int,T,L,L,2,2)
        data["anyon_hist"] = zeros(Bool,T,L,L)  

        # create an example fixed error pattern 
        state = falses(L,L,2)
        hL = div(L,2)
        dx = div(L,5)
        state[hL-dx:hL+dx,hL,1] .= true
        state[hL+dx+1,hL:hL+dx,2] .= true
        anyons = get_synds(state)
        fields = zeros(Int,L,L,2,2)
        new_fields = zeros(Int,L,L,2,2)
        correction = falses(L,L,2)

        ### do a search for configurations that take a long time to decay ### 
        found_example = false 
        examples = 0 
        max_examples = 1
        while ~found_example && examples < max_examples
            examples += 1 
            fields .*= 0
            data["field_hist"] .= 0 
            data["hist"] .= false
            data["anyon_hist"] .= false
            if examples > -1
                state,anyons = init_2d(p,L,"rand",bconds)
            end 
        
            for t in 1:T
                data["hist"][t,:,:,:] .= state
                data["field_hist"][t,:,:,:,:] .= fields
                data["anyon_hist"][t,:,:] .= anyons
                update_state!(state,correction,anyons,fields,new_fields,r,0,synch,bconds)
                if ~any(anyons) # if we are anyon-free 
                    break 
                end
            end 
            if examples % 100 == 0 println(examples) end 
            if any(anyons) 
                found_example = true 
            end
        end 

    ### offline decoding ### 
    elseif mode == "erode" 
        println("doing offline decoding at L = $L...")
        println("boundary conditions = $bconds")
        data["erode_frac"] = zeros(nps)
        data["erode_times"] = zeros(nps)
        data["erode_stats"] = zeros(nps,3) # mean, std, max of erosion times
        C = 2 
        # data["anyon_density"] = zeros(nps,C*maximum(Ls)) # average anyon density at time t (not super useful; currently not calculated)
        for (pind,thisp) in enumerate(ps) 
            thisL = Ls[pind]
            max_erode_time = C * thisL^2
            println("L = $thisL, p = $thisp")
            longest_erosion = 0 
            fields = zeros(Int,thisL,thisL,2,2); new_fields = zeros(Int,thisL,thisL,2,2); correction = falses(thisL,thisL,2)
            @showprogress dt=1 for _ in 1:samps[pind] 
                # create suitably random initial state # 
                state,anyons = init_2d(thisp,thisL,"rand",bconds) # second returned argument is charges, which we don't need to make use of 
                t = 0
                while t < max_erode_time  # no absorbing steady states since particles move stochastically
                    update_state!(state,correction,anyons,fields,new_fields,r,0,synch,bconds)
                    t += 1 
                    if any(anyons) # if we are not in a logical state 
                        # data["anyon_density"][pind,t] += sum([state[i] ⊻ state[mod1(i+1,thisL)] for i in 1:thisL]) / thisL  
                        continue
                    else # if the state is charge-free, erosion is done 
                        break 
                    end 
                end 
                if t == max_erode_time 
                    println("max erosion time reached!") 
                end 
                data["erode_frac"][pind] += detect_logical_error(state) / samps[pind]

                if t > longest_erosion longest_erosion = t end

                data["erode_times"][pind] += t / samps[pind]
                data["erode_stats"][pind,1] += t / samps[pind]
                data["erode_stats"][pind,2] += t^2 / samps[pind]
                # data["anyon_density"][pind,:] ./= samps[pind] 
            end 
            data["erode_stats"][pind,2] = sqrt(data["erode_stats"][pind,2] - data["erode_stats"][pind,1]^2 + 1e-10)
            data["erode_stats"][pind,3] = longest_erosion
            println("⟨t⟩ = ",data["erode_stats"][pind,1])
            println("σ(t) = ",data["erode_stats"][pind,2])
            println("max(t) = ",data["erode_stats"][pind,3])
        end 

    else ## stuff requiring monte carlo averages

        ### initialize various things ### 
        nsteps = vary_L ? length(Ls) : length(ps)
        scalar_quantities = ["Ft" "binds" "chis" "Ms" "trels"]
        for key in scalar_quantities data[key] = zeros(nsteps) end
        data["trel_stats"] = zeros(nsteps,3) # mean, std, max of relaxation times

        ### compute things requiring monte carlo averages ### 
        function compute(p,samps,L,T)
            """ 
            computes monte carlo averages of various quantities for a given value of (p,q,L)
            returns: 
                mode = trel: relaxation time
                mode = Ft: decoding fidelity, time-dependent magnetization
            uses synch, eta as global parameters (only p, q, and L may vary)
            """

            mc_keys = ["trels" "Ft" "trel_stats"] # quantities to compute
            mc_data = Dict{String, Any}(key => 0 for key in mc_keys) # store the results of the monte carlo averages
            # mc_data["Mt"] = zeros(T) # time-dependent magnetization

            correction = falses(L,L,2); new_fields = zeros(Int,L,L,2,2) 
            state = falses(L,L,2); fields = zeros(Int,L,L,2,2); anyons = falses(L,L) 

            if mode == "trel"
                mc_data["trel_stats"] = zeros(3) # mean, std, max of relaxation times
                maxT = 50000000
                logical_mag = 1 # begin the system in the logical state aligned against the bias of the noise 
                max_trel = 0 
                decode_attempt_period = 5 
                fields_to_decode = zeros(Int,L,L,2,2)
                @showprogress dt=1 desc="sampling..." for samp in 1:samps 
                    state .= false; anyons .= false; fields .*= 0 
                    t = 1 
                    while t < maxT
                        update_state!(state,correction,anyons,fields,new_fields,r,p,synch,bconds)
                        if t%decode_attempt_period == 0 
                            state_to_decode = copy(state); anyons_to_decode = copy(anyons)
                            fields_to_decode .*= 0
                            for _ in 1:T 
                                update_state!(state_to_decode,correction,anyons_to_decode,fields_to_decode,new_fields,r,0,synch,bconds)
                                if ~any(anyons_to_decode) 
                                    break 
                                end
                            end
                            if ~detect_logical_error(state_to_decode) break end 
                        end 
                        t += 1 
                    end
                    if t == maxT 
                        println("maxT reached! (sample = $samp)")
                    end
                    mc_data["trels"] += t / samps 
                    if t > max_trel max_trel = t end
                    mc_data["trel_stats"][1] += t / samps
                    mc_data["trel_stats"][2] += t^2 / samps
                end 
                mc_data["trel_stats"][2] = sqrt(mc_data["trel_stats"][2] - mc_data["trel_stats"][1]^2)
                mc_data["trel_stats"][3] = max_trel
            end 

            if mode == "Ft"
                println("T = $T")
                @showprogress dt=1 desc="sampling..." for s in 1:samps 
                    state .= true  
                    fields .*= 0 
                    this_sum = sum(state)
                    
                    for _ in 1:T 
                        update_state!(state,correction,fields,new_fields,r,p,synch)
                        this_sum = sum(state)
                    end 
                    mc_data["Ft"] += this_sum > L/2 ? 1/samps : 0 # decoding fidelity at end time
                end 
            end

            return mc_data
        end 

        for i in 1:nsteps 

            thisp = p; thisL = L; thissamps = samps[i]
            thisT = Ts[1]
            if vary_L   
                thisL = Ls[i]; thisT = Ts[i]
                println("L = $thisL | samps : $thissamps")
            else
                thisp = ps[i]
                println("p = $thisp | samps : $thissamps")
            end
            this_mc_data = compute(thisp,thissamps,thisL,thisT)
            for key in keys(this_mc_data)
                if key == "trel_stats"
                    data[key][i,:] = this_mc_data[key]
                    println("trel stats: ",this_mc_data[key])
                else
                    data[key][i] = this_mc_data[key]
                end 
            end
        end 
    end 

    # write to file 
    sadj = ~synch ? "_asynch" : ""    
    fout = "tmp_data/2dloc_$(mode)_"*(vary_L ? "p$p" : "L$L")*"$sadj$out_adj.jld2"

    println("writing to file: $fout")
    f = jldopen(fout,"w")
    for key in keys(data)
        write(f,key,data[key])
    end 
    for key in keys(params)
        write(f,key,params[key])
    end 
    close(f)    
    if vary_L
        alert("finished | p = $p; L = $(Ls[1]) → $(Ls[end])")
    else 
        alert("finished | L = $L; p = $(ps[1]) → $(ps[end])")
    end 
    println("finished at time $(Dates.now())")
end 

main()