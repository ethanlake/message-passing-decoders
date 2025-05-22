include("common_functions.jl")

"""
simulation of the 1D local message passing decoder
"""

function update_fields!(fields,new_fields,anyons)
    L,_ = size(fields)
    hL = div(L,2)
    ind(i) = mod1(i,L)

    for i in 1:L 
        ip1 = ind(i+1); im1 = ind(i-1)
        # left-moving fields
        new_fields[i,1] = anyons[ip1] ? 1 : (fields[ip1,1] == 0 ? 0 : fields[ip1,1] + 1) # left propagating fields
        # right-moving fields
        new_fields[i,2] = anyons[im1] ? 1 : (fields[im1,2] == 0 ? 0 : fields[im1,2] + 1) # right propagating fields
    end 
    new_fields[findall(new_fields .> L)] .= 0 
    fields .= new_fields
end 

function update_state!(state,correction,fields,new_fields,r,p,synch)
    """
    state, correction: L of Bools 
    fields, new_fields: Lx2 of Ints  
    r: ratio of field updates to spin updates
    p: error probability 
    synch: if true, updates are done synchronously 
    """

    L,_ = size(fields)
    hL = L ÷ 2
    ind(i) = mod1(i,L)

    dws = state .⊻ circshift(state,-1) # dws[i] = (domain wall between i and i+1)

    if synch 
        new_fields .*= 0 
        for _ in 1:r  
            update_fields!(fields,new_fields,dws)
        end 

        # update (all) the anyon positions 
        correction .= false 
        for i in 1:L  
            if dws[i] # domain wall between sites i,i+1 
                if fields[i,1] != 0 && fields[i,2] != 0 
                    force = -sign(fields[i,1] - fields[i,2]) # if fields[i,1] > fields[i,2] and both nonzero, anyon to the right is further than anyon to the left; thus move left if force < 0. 
                    # if rand() < .9 # some stochastic part when forces are balanced to break out of limit cycles at small system sizes 
                        if force > 0 
                            correction[ind(i+1)] = true
                        elseif force < 0 
                            correction[i] = true
                        end 
                    # else 
                    #     dir = rand([0,1])
                    #     correction[ind(i+dir)] = true
                    # end     
                elseif fields[i,1] != 0 # only left-moving messages
                    correction[ind(i+1)] = true
                elseif fields[i,2] != 0 # only right-moving messages
                    correction[i] = true
                end 
            end 
        end
        state .⊻= correction

        # add noise 
        for i in 1:L 
            if rand() < p 
                state[i] ⊻= true 
            end 
        end

    ### asychronous updates ### 
    else 
        for _ in 1:round(L/2) # poocontinuoustime
            i = rand(1:L); ip1 = ind(i+1); im1 = ind(i-1)

            ### update fields ### 
            fields[i,1] = dws[ip1] ? 1 : (fields[ip1,1] == 0 ? 0 : fields[ip1,1] + 1) # left propagating fields
            fields[i,2] = dws[im1] ? 1 : (fields[im1,2] == 0 ? 0 : fields[im1,2] + 1) # right propagating fields

            if fields[i,1] > L fields[i,1] = 0 end
            if fields[i,2] > L fields[i,2] = 0 end 

            ### update anyon positions ### 
            if dws[i] # domain wall between sites i,i+1 
                if fields[i,1] != 0 && fields[i,2] != 0 
                    force = -sign(fields[i,1] - fields[i,2]) # if fields[i,1] > fields[i,2] and both nonzero, anyon to the right is further than anyon to the left; thus move left if force < 0. 
                    if force > 0 # move dw right
                        state[ind(i+1)] ⊻= true
                        dws[i] ⊻= true; dws[ip1] ⊻= true
                    elseif force < 0 # move dw left 
                        state[i] ⊻= true
                        dws[i] ⊻= true; dws[im1] ⊻= true
                    end     
                elseif fields[i,1] != 0 # only left-moving messages
                    state[ip1] ⊻= true
                    dws[i] ⊻= true; dws[ip1] ⊻= true
                elseif fields[i,2] != 0 # only right-moving messages
                    state[i] ⊻= true 
                    dws[i] ⊻= true; dws[im1] ⊻= true
                end 
            end 

        end 
    end 
end 

function parameter_repository(mode,L,p,vary_L)
    pmin = 0; pmax = 1; nps = 1; nLs = 1 
    samps = 1; samps_vec = [samps]; ps = [p]; Ls = [L]; Ts = [1]

    if mode == "trel"
        ### fix system size, vary p ###
        if ~vary_L 
            nps = 10
            samps = round(Int,20000 / L)
            samps_vec = [samps for _ in 1:nps]
            
            pmin = .07; pmax = .23
            # if L > 16 pmin = .06 end 
            # if L > 32 pmin = .07 end 
            ps = [10^(x) for x in LinRange(log10(pmin),log10(pmax),nps)] 

        ### vary system size, fix p ###
        else 
            Lmin = 8; Lmax = 2*1024
            # if p < .15 Lmax = 512 end 
            # if p ≥ .135 Lmax = 2048 end
            nps = 15
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nps)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,4500000/el) for el in Ls] # poosamps 
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
            samps_vec = [round(Int,800000/el) for el in Ls] # poosamps 
            Ts = [el for _ in Ls]
        end 
    end 

    if mode == "erode"
        
        if ~vary_L 
            nps = 15
            samps_vec = [round(Int,5000000 / L) for _ in 1:nps]
            pmin = .3; pmax = .5 
            ps = [p for p in LinRange(pmin,pmax,nps)]

        else # varying L 
            # for looking at the erosion time: can go to very large L 
            Lmin = 800; Lmax = 48*2056 
            # Lmax = 16 * 2056
            nps = 9
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nps)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,100620000/el) for el in Ls] # professional 
            samps_vec = [round(Int,10620000/el) for el in Ls] # poosamps 

            # for getting the logical error rate: use smaller system sizes 
            Lmin = 16; Lmax = 256 
            if p < .45 
                Lmax = 175
            end 
            if p < .35 
                Lmax = 100 
            end 
            nps = 6 
            Ls = [round(Int,el) for el in LinRange(Lmin,Lmax,nps)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,525000000/el) for el in Ls] # poosamps 
            samps_vec = [round(Int,85000000/el) for el in Ls] # poosamps 

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

function adversarial_state_gen(n,k,l,seed)
    """
    generates states using an l-fold iteration of the recursion relation 
    0 -> 0^n 
    1 -> 1^k 0^{n-k}
    on input state seed
    """
    rules = Dict(
        0 => fill(0, n),
        1 => vcat(fill(1, k), fill(0, n - k))
    )
    
    # Initialize with the starting symbol
    current = [seed ? 1 : 0]
    
    # Apply the substitution rule l times
    for _ in 1:l
        next = Int[]
        for symbol in current
            append!(next, rules[symbol])
        end
        current = next
    end
    
    return Bool.(current)
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
    L = 1028
    L = 48
    L = 16
    L = 200
    L = 96
    # pool 
    L = 8
    L = 16 
    L = 32
    L = 256 
    L = 64
    L = 128 
    L = 256
    L = 512
    p = .35 # poop (used if varying L)
    vary_L = ~true # if true, vary system size; if false, use fixed system size and vary p # poovaryL 

    r = 3 # number of field updates per spin update 
    synch = ~true 

    out_adj = ""

    params = parameter_repository(mode,L,p,vary_L)
    Ts = params["Ts"]; samps = params["samps"]; ps = params["ps"]; nps = params["nps"]; p = params["p"];  Ls = params["Ls"]
    params["synch"] = synch; params["r"] = r

    data_keys = ["Ft" "Mt"] ∪ ["hist" "field_hist"] ∪ ["trels" "trel_stats"] ∪ ["erode_times" "erode_stats" "anyon_density"] 
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

    state = falses(L); correction = falses(L); fields = zeros(Int,L,2); new_fields = zeros(Int,L,2)

    ### write history of evolution ### 
    if mode == "hist" 
        T = round(Int,1L) 
        println("running for time T = $T")
        data["hist"] = zeros(Bool,T,L)
        data["field_hist"] = zeros(Int,T,L,2)
        state[20:40] .= true
        # state = rand(Bool,L)
        
        ### test out some adversarial situations ### 
        # n = 6; k = n-1; l = 6
        # r = 3
        # L = n^l 
        # state = adversarial_state_gen(n,k,l,true)
        # println("initial state = $state")
        # fields = zeros(L,2); new_fields = zeros(Int,L,2); correction = falses(L)
        # terode = 0 
        # while ~all(state) && any(state)
        #     terode += 1 
        #     update_state!(state,correction,fields,new_fields,r,0,synch)
        # end 
        # println("final log bit = $(Int(state[1]))")
        # println("erosion time / L = $(terode/L)")

        # return 

        ### do a search for configurations that take a long time to decay ### 
        found_example = false 
        examples = 0 
        max_examples = 1 
        while ~found_example && examples < max_examples
            examples += 1 
            fields = zeros(L,2); new_fields = zeros(Int,L,2)
            if examples > 0 state = init(.5,L,"rand") end 
        
            for t in 1:T
                data["hist"][t,:] .= state
                data["field_hist"][t,:,:] .= fields
                update_state!(state,correction,fields,new_fields,r,0,synch)
                # if all(state) || ~any(state) # if we are anyon-free 
                #     break 
                # end
            end 
            if examples % 100 == 0 println(examples) end 
            if any(state) && ~all(state) 
                found_example = true 
            end
        end 

    ### offline decoding ### 
    elseif mode == "erode" 
        println("doing offline decoding at L = $L...")
        data["erode_frac"] = zeros(nps)
        data["erode_times"] = zeros(nps)
        data["erode_stats"] = zeros(nps,3) # mean, std, max of erosion times
        C = 2 
        data["anyon_density"] = zeros(nps,C*maximum(Ls)) # average anyon density at time t 
        for (pind,thisp) in enumerate(ps) 
            thisL = Ls[pind]
            fields = zeros(Int,thisL,2); new_fields = zeros(Int,thisL,2)
            correction = falses(thisL); state = falses(thisL)
            max_erode_time = thisL^2 
            println("L = $thisL, p = $thisp")
            longest_erosion = 0 
            @showprogress dt=1 for _ in 1:samps[pind] 
                # create suitably random initial state # 
                state = init(thisp,thisL,"rand") # second returned argument is charges, which we don't need to make use of 
                fields .*= 0 # begin with trivial fields and forces 
                init_logical = sum(state) < thisL/2 ? 0 : 1
                t = 0
                while t < max_erode_time  # no absorbing steady states since particles move stochastically
                    update_state!(state,correction,fields,new_fields,r,0,synch)
                    t += 1 
                    if any(state) && !all(state) # if we are not in a logical state 
                        # data["anyon_density"][pind,t] += sum([state[i] ⊻ state[mod1(i+1,thisL)] for i in 1:thisL]) / thisL  
                        continue
                    else # if the state is charge-free, erosion is done 
                        break 
                    end 
                end 
                if t == max_erode_time 
                    println("max erosion time reached!") 
                end 
                data["erode_frac"][pind] += (state[1] == init_logical ? 1 : 0) / samps[pind] 

                if t > longest_erosion longest_erosion = t end

                data["erode_times"][pind] += t / samps[pind]
                data["erode_stats"][pind,1] += t / samps[pind]
                data["erode_stats"][pind,2] += t^2 / samps[pind]
                data["anyon_density"][pind,:] ./= samps[pind] 
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

            correction = falses(L); new_fields = zeros(Int,L,2) 
            state = falses(L); fields = zeros(Int,L,2) 

            if mode == "trel"
                mc_data["trel_stats"] = zeros(3) # mean, std, max of relaxation times
                maxT = 50000000
                logical_mag = 1 # begin the system in the logical state aligned against the bias of the noise 
                max_trel = 0 
                @showprogress dt=1 desc="sampling..." for samp in 1:samps 
                    state .= true 
                    fields .*= 0 
                    t = 1 
                    while t < maxT
                        update_state!(state,correction,fields,new_fields,r,p,synch)
                        if 2*(sum(state)/L - .5) ≤ 0 # logical error 
                            break 
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
            thisT = Ts[i]
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
                else
                    data[key][i] = this_mc_data[key]
                end 
            end
        end 
    end 

    # write to file 
    sadj = ~synch ? "_asynch" : ""    
    fout = "tmp_data/1dloc_$(mode)_"*(vary_L ? "p$p" : "L$L")*"$sadj$out_adj.jld2"

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
    # print time that script finished 
    println("finished at time $(Dates.now())")
end 

main()