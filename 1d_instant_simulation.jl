include("common_functions.jl")

"""
automaton where domains produce fields with power law force interactions of the form 1/x^alpha. 
can set dw motion to either have a bias like 1/x^alpha, or to always move in the direction where the force is largest. 
domain walls come with charges, which can either be assigned randomly (as they would be for charged anyons in 2D) or in an alternating fashion (as they would be for dws in 1D)
"""

function update_state!(state,charges,charged,p,q,synch,alpha,gamma,interactions,correction_field,correction,charge_correction;alter=false)
    """
    state: Bool of spins 
    charges: charges[i] is the charge of the domain wall between sites i,i+1 
    charged: wether or not to do charged decoding (if not, charges is vector of 0s and 1s, if yes, is a vector of -1s, 0s, and 1s that sums to 0)
    p: error probability 
    q: measurement error probability 
    alpha: power law exponent for force interaction between domain walls
    gamma: strength of the bias (relative to diffusion). domain wall moves left with probability (1+v)/2, where v = gamma F / (gamma + |F|). if gamma = 0, then the domain wall moves in the direction of the force (indepdent of the force's magnitude)
    alter: if true, charged decoding is done using charges of alternating sign (as dws would be)
    """
    L = size(state)[1]
    hL = L ÷ 2
    ind(i) = mod1(i,L)
    qcharged = charged ? 1 : -1 # uncharged case is attractive; this sign is needed to simplify things 
    if charged @assert q == 0 "measurement errors not implemented yet"; @assert sum(charges) == 0 "not charge neutral!" end 

    @assert ~charged "something is fucked up about charged updates... check history with bias towards false"

    if synch                 
        @views dws = findall(charges .!= 0)
            
        # reset the correction fields 
        fill!(correction_field, 0.0); fill!(charge_correction, 0)

        for dw in dws
            qi = charges[dw] # calculate the force field produced by this domain wall on all sites except the one it occupies 
            for j in 1:hL 
            # set up the force fields produced by this domain wall
                if charges[ind(dw+j)] != 0 
                    correction_field[ind(dw+j)] += qi * qcharged * interactions[dw,ind(dw+j)] #circle_distance(dw,dw+j,L)^alpha # this dw wants to make those to its right move left if they have opposite charges or if we are doing uncharged decoding
                end
                if charges[ind(dw-j)] != 0 
                    correction_field[ind(dw-j)] -= qi * qcharged * interactions[dw,ind(dw-j)] #circle_distance(dw,dw-j,L)^alpha # this dw wants to make those to its left move to the right
                end
                # correction_field[ind(dw+j)] += qi * qcharged * interactions[dw,ind(dw+j)] #circle_distance(dw,dw+j,L)^alpha # this dw wants to make those to its right move left if they have opposite charges or if we are doing uncharged decoding
                # correction_field[ind(dw-j)] -= qi * qcharged * interactions[dw,ind(dw-j)] #circle_distance(dw,dw-j,L)^alpha # this dw wants to make those to its left move to the right
            end
        end 

        # now move each domain wall in the direction of the force it feels 
        rands = rand(length(dws))
        if charged 
            for i in randperm(length(dws)) 
                dw = dws[i]
                qi = charges[dw]
                if qi != 0 #&& rand() < .9 # can add laziness to remove unwanted limit cycles but perhaps not a big deal?  
                    F = qi * correction_field[dw] 
                    v = F / (abs(F) + gamma)
                    if rands[i] ≤ (1+v)/2 # if dw wants to move right 
                        if charges[ind(dw+1)] != charges[dw]
                            state[ind(dw+1)] ⊻= true
                            charges[dw] -= qi 
                            charges[ind(dw+1)] += qi
                        end 
                    else # if dw wants to move left 
                        if charges[ind(dw-1)] != charges[dw]
                            state[dw] ⊻= true 
                            charges[dw] -= qi
                            charges[ind(dw-1)] += qi
                        end
                    end
                end 
            end 
        else 
            fill!(correction, false)
            for i in 1:length(dws) 
                dw = dws[i]
                qi = charges[dw]
                if qi != 0 && rand() < .9 # just adding laziness to remove unwanted limit cycles but perhaps not a big deal?
                    F = qi * correction_field[dw] 
                    v = F / (abs(F) + gamma)
                    if rands[i] ≤ (1+v)/2 # if dw wants to move right 
                        correction[ind(dw+1)] = true
                    else # if dw wants to move left 
                        correction[dw] = true 
                    end
                end 
            end 
            state .⊻= correction
        end 

        # now do the errors -- happens in series for different sites since have to sort out the charge updates
        if charged 
            rands = rand(L)
            for i in randperm(L) # more efficient than shuffle(1:L) which creates a temp array 
                if rands[i] < p 
                    if alter    
                        state[i] ⊻= true # spin flip always allowed since charges alternate; charges will be determined at the end 
                    else 
                        qn = rand([-1,1]) # try to add anyon of charge charge at site i, and charge -charge at site i+1. if this is not possible, try to do the same with the opposite sign of charge. 
                        if charges[i] != qn  && charges[ind(i-1)] != -qn 
                            state[i] ⊻= true 
                            charges[i] += qn ; charges[ind(i-1)] -= qn 
                        elseif -charges[i] != qn && -charges[ind(i-1)] != -qn 
                            state[i] ⊻= true 
                            charges[i] -= qn ; charges[ind(i-1)] += qn 
                        end
                    end 
                end
            end      
            if alter # determine the charges just from the alternating signs of the domain walls 
                for i in 1:L
                    if state[i] != state[ind(i+1)]
                        charges[i] = state[i] ? 1 : -1 # ensures the charges are always alternating in sign 
                    else 
                        charges[i] = 0 
                    end 
                end 
            end
        else 
            rands = rand(L)
            for i in 1:L 
                if rands[i] < p 
                    state[i] ⊻= true 
                end
            end
            for i in 1:L # note: doing this as charges = [state != state for i in 1:L] would NOT modify charges after the function call finished 
                charges[i] = state[i] != state[ind(i+1)] ? 1 : 0  
            end 
        end
    else 
        println("asynch not implemented yet")
    end     
end 

function parameter_repository(mode,L,p,vary_L,measurement_error_strength,synch,alpha,gamma,charged)
    q = measurement_error_strength * p 
    pmin = 0; pmax = 1; nps = 1; nLs = 1 
    samps = 1; ps = [p]; qs = [measurement_error_strength * p]; Ls = [L]
    T = 1; Ts = [1]
    samps_vec = [samps] # number of samples used for each parameter value (p,L)

    pc = 0

    if mode == "trel"

        ### fix system size, vary p ###
        if ~vary_L 
            nps = 10

            samps = round(Int,45000 / (L))
            T = 1

            if charged 
                if alpha == .25  # p = .25 seems to be slightly below threshold 
                    if L < 48
                        pmin = .16; pmax = .45 # .15 to .23 at T = L if doing alternating domains 
                    elseif L < 64 
                        pmin = .21; pmax = .45 # .15 to .23 at T = L if doing alternating domains 
                    else 
                        pmin = .25; pmax = .45
                    end 
                elseif alpha == .75
                    pmin = .01; pmax = .1 
                elseif alpha == 1.5 
                    pmin = .01; pmax = .1 
                elseif alpha == 2.5 
                    pmin = .01; pmax = .1 
                end 
            else 
                if alpha == .25 # by p ~ .05 trel actually decreases with L. clearly lacks a threshold
                    pmin = .015; pmax = .08 # .15 to .23 at T = L if doing alternating domains 
                    pmin = .015; pmax = .25 # .15 to .23 at T = L if doing alternating domains 
                elseif alpha == .5 # trel vs p indicates weak transition around .03
                    pmin = .02; pmax = .12
                elseif alpha == .75 # again pc ~ .05 from trel vs p, but seems very weak --- probably no threshold but need to check trel vs L at p ~ .05 to be sure  
                    pmin = .035; pmax = .135 
                    # if L > 100 pmin = .04 end 
                elseif alpha == 1 # pc ~ .05 maybe from trel vs p 
                    pmin = .04; pmax = .1
                elseif alpha == 1.5 # pc ~ .05 / .06 from trel vs p, but very weak -- probably no threshold 
                    if L == 32 
                        pmin = .035; pmax = .125
                    end 
                    if L == 64 
                        pmin = .039; pmax = .087
                    end 
                    if L == 128 
                        pmin = .038; pmax = .06
                        nps = 6; samps = round(Int,5000 / (L))
                    end
                elseif alpha == 2 
                    pmin = .06; pmax = .16
                elseif alpha == 2.5 
                    pmin = .01; pmax = .1 
                end 
            end 

            samps_vec = [samps for i in 1:nps]

        ### vary system size, fix p ###
        else 
            Lmin = 32; Lmax = 512
            nels = 10
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nels)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,800000/el) for el in Ls] # poosamps 
            Ts = [1 for _ in Ls]

        end
    end 

    if mode == "Ft"
        nps = 11
        samps_vec = [round(Int,200000 / (L)) for _ in 1:nps]

        if ~vary_L  # vary p 
            Ts = [L]; T = L
            if charged 
                if alpha == .25 
                    pmin = .35; pmax = .8 # .15 to .23 at T = L if doing alternating domains 
                elseif alpha == .75
                    pmin = .01; pmax = .1 
                elseif alpha == 1.5 
                    pmin = .01; pmax = .1 
                elseif alpha == 2.5 
                    pmin = .01; pmax = .1 
                end 
            else 
                if alpha == .25 # clearly (?) lacks a threshold
                    pmin = .01; pmax = .1 # .15 to .23 at T = L if doing alternating domains 
                elseif alpha == .75
                    pmin = .01; pmax = .1 
                elseif alpha == 1.5 # 
                    pmin = .01; pmax = .1 
                elseif alpha == 2.5 
                    pmin = .01; pmax = .1 
                end 
            end 

        else # vary L  
            if alpha == .25 
                pmin = .35; pmax = .8 # .15 to .23 at T = L if doing alternating domains 
            elseif alpha == .75
                pmin = .01; pmax = .1 
            elseif alpha == 1.5 
                pmin = .01; pmax = .1 
            elseif alpha == 2.5 
                pmin = .01; pmax = .1 
            end 

        end 

    end 

    if mode == "erode"
        nps = 10
        samps_vec = [round(Int,1000000 / (L)) for _ in 1:nps]
        pmin = .01; pmax = .5 

        if ~vary_L 
            if charged 
                if alpha == .25 
                    pmin = .05; pmax = .5 # .15 to .23 at T = L if doing alternating domains 
                elseif alpha == .75
                    pmin = .01; pmax = .1 
                elseif alpha == .5 
                    pmin = .01; pmax = .5 
                    # pmin = .125; pmax = .35
                elseif alpha == 1 
                    pmin = .01; pmax = .1
                elseif alpha == 1.5 
                    pmin = .01; pmax = .5 
                elseif alpha == 2.5 
                    pmin = .01; pmax = .1 
                end 
            else # uncharged
                
                if alpha == .01 
                    pmin = .01; pmax = .5 
                elseif alpha == .25 
                    pmin = .01; pmax = .5 # .15 to .23 at T = L if doing alternating domains 
                elseif alpha == .5 
                    pmin = .01; pmax = .5
                elseif alpha == .75
                    pmin = .25; pmax = .5
                elseif alpha == 1. 
                    pmin = .01; pmax = .5
                elseif alpha == 1.5 
                    pmin = .01; pmax = .4 
                elseif alpha == 2.5 
                    pmin = .01; pmax = .1 
                end 
            end 
        
        else # varying L 
            Lmin = 16; Lmax = 2*2056; nps = 14
            Ls = [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nps)]
            Ls .+= Ls .% 2 # make sure Ls are even to avoid even / odd effect 
            reverse!(Ls) # do the long ones first so we have an idea of how long things will take
            samps_vec = [round(Int,2200000/el) for el in Ls] # poosamps 
            Ts = [1 for _ in Ls]
        end 
    end 

    if mode ∈ ["erode_times" "hist"]
        ps = [pmin]; nps = 1  
    else
        if vary_L 
            ps = [p for _ in 1:length(Ls)]
            nps = length(Ls)
        else 
            ps = [p for p in LinRange(pmin,pmax,nps)]
            if mode ∈ ["trel"] # space out evenly on a log scale instead 
                ps = [10^(x) for x in LinRange(log10(pmin),log10(pmax),nps)] 
            end 
        end 
    end 
    qs = measurement_error_strength .* ps  

    params = Dict{String, Any}()
    params["samps"] = samps_vec; params["T"] = T; params["Ts"] = Ts; params["mode"] = mode; params["L"] = L; params["nps"] = nps; params["ps"] = ps; params["qs"] = qs; params["synch"] = synch; params["pc"] = pc; params["Ls"] = Ls; params["p"] = p; params["q"] = q; params["vary_L"] = vary_L; params["alpha"] = alpha; params["gamma"] = gamma; params["charged"] = charged

    return params 
end 

function main()

    """
    supported modes: 
    * "stats": get statistics in the long time steady states 
    * "trel": compute relaxation time
    * "Ft": get decoding fidelity and time-dependent magnetization  
    * "hist": history of particular evolution
    * "erode": check if initial errors get eroded under ideal evolution 
    * "erode_times": computes erosion time of initial domain walls of different sizes 
    """

    mode = "trel" # "stats" "trel" "Ft" "hist" "erode" "erode_times" # poomode 
    # pool # 48 64 96 128 192 256 384 512 (used if varying p)
    L = 48
    L = 256
    L = 1028
    L = 128
    L = 512
    L = 8 
    L = 16
    L = 32 
    L = 64
    L = 128 
    
    p = .01 # poop (used if varying L)
    vary_L = ~true # if true, vary system size; if false, use fixed system size and vary p # poovaryL 
    ind(i) = mod1(i,L)

    alpha = 2. # power law exponent  # pooalpha 
    charged = ~true # poocharge 
    alter = false # if charged, makes domain walls alternate in signs (as they would for dws in 1D); if uncharged, makes domain walls have random charges (as they would for charged anyons in 2D) # pooalter 
    gamma = 0 # bias of noise is force / (gamma + |force|) --- becomes sgn(F) in the limit gamma = 0 # poogamma
    synch = true # whether or not to synchronize the updates # poosynch 

    out_adj = "" # poooutadj
    measurement_error_strength = 0 # syndrome measurements fail with probability measurement_error_strength * p # pooq 
    q = p * measurement_error_strength

    params = parameter_repository(mode,L,p,vary_L,measurement_error_strength,synch,alpha,gamma,charged)
    Ts = params["Ts"]; T = params["T"]; samps = params["samps"]; ps = params["ps"]; qs = params["qs"]; nps = params["nps"]; p = params["p"];  Ls = params["Ls"]

    data_keys = ["erode_sizes" "erode_frac" "mags"] ∪ ["Ft" "Mt"] ∪ ["Ms" "chis" "binds"] ∪ ["hist" "field_hist" "charge_hist"] ∪ ["trels"] ∪ ["erode_times"]   # quantities to store
    data = Dict{String, Any}(key => 0 for key in data_keys)

    println("details of simulation: ")
    if vary_L 
        println("(p,q) = ($p,$q)")
        println("Ls = $Ls")
        # println("Ts = $Ts")
    else 
        println("system size: $L")
        println("ps = $ps")
        println("measurement error strength = $measurement_error_strength")
    end 
    println("synch = $synch")
    println("charged = $charged")
    println("mode = $mode")
    println("α = $alpha")   
    println("γ = $gamma")

    # calculate distance functions to be used when updating state 
    interactions = zeros(L,L)
    for i in 1:L 
        for j in 1:L 
            if i != j
                interactions[i,j] = 1/circle_distance(i,j,L)^alpha 
                # interactions[i,j] = exp(-5.0 * circle_distance(i,j,L)) # to test nearest-dw interactions
            end 
        end 
    end
    # define correction fields which will be passed as arguments to update function 
    correction_field = zeros(L)
    correction = zeros(Bool,L)
    charge_correction = zeros(Int,L)

    ### write history of evolution ### 
    if mode == "hist" # poohist 
        println("running history at size $L with noise (p,q) = ($(ps[1]),$(qs[1]))...")
        Ts[1] = round(Int,L/2)
        println("running for time T = $(Ts[1])")
        data["hist"] = zeros(Bool,Ts[1],L)
        data["field_hist"] = zeros(Ts[1],L)
        data["charge_hist"] = zeros(Int,Ts[1],L)

        failure = false
        maxcounts = 1
        counts = 0
        state = trues(L); charges = zeros(Int,L)
        while ~failure && counts < maxcounts 
            state,charges = charged_init(.6,L,charged,alter,"rand")
            init_log = sum(state) < L/2 ? 0 : 1 
            # println("init_log = $init_log")
            for t in 1:Ts[1]
                data["hist"][t,:] .= state
                data["charge_hist"][t,:] .= charges
                update_state!(state,charges,charged,.0,qs[1],synch,alpha,gamma,interactions,correction_field,correction,charge_correction;alter=alter)
                data["field_hist"][t,:] .= correction_field
            end 
            if (sum(state) < L/2 ? 0 : 1) != init_log
                failure = true 
            end
            counts += 1 
        end 
    
    elseif mode == "erode"   # random initial states with errors determined by p 
        println("doing error-free erosion at L = $L...")
        data["erode_frac"] = zeros(nps)
        data["erode_times"] = zeros(nps)
        nsteps = vary_L ? length(Ls) : length(ps); @assert nps == nsteps 

        for (pind,thisp) in enumerate(ps) 
            thisL = vary_L ? Ls[pind] : L
            println("p = $thisp, L = $thisL")
            if vary_L # precompute the interaction kernel for this system size 
                interactions = zeros(thisL,thisL)
                for i in 1:thisL 
                    for j in 1:thisL 
                        if i != j
                            interactions[i,j] = 1/circle_distance(i,j,thisL)^alpha 
                            # interactions[i,j] = exp(-5.0 * circle_distance(i,j,thisL)) # to test nearest-dw interactions 
                        end 
                    end 
                end
                # define correction fields which will be passed as arguments to update function 
                correction_field = zeros(thisL)
                correction = zeros(Bool,thisL)
                charge_correction = zeros(Int,thisL)
            end             
            
            @showprogress dt=1 for _ in 1:samps[pind] 

                # create suitably random initial state # 
                state,charges = charged_init(thisp,thisL,charged,alter,"rand")
                init_logical = sum(state) < thisL/2 ? 0 : 1
                # println("mag = ",2*(sum(state)/L - .5))
                t = 0 
                # println("b: ",sum(abs.(charges)))
                tmax = 50thisL^2
                while t < tmax # finite time allowing for possibility of other absorbing states 
                    update_state!(state,charges,charged,0,0,synch,alpha,gamma,interactions,correction_field,correction,charge_correction;alter=alter)
                    # println("a: ",sum(abs.(charges)))

                    t += 1 
                    if any(charges .!= 0) #
                        continue
                    else # if the state is charge-free, erosion is done 
                        break 
                    end 
                end 
                if t == tmax println("max time reached!") end 
                data["erode_frac"][pind] += (state[1] == init_logical ? 1 : 0) / samps[pind] 
                data["erode_times"][pind] += t / samps[pind]
            end 
        end 

    else ## stuff requiring monte carlo averages (other than erode_times)

        ### initialize various things ### 
        nsteps = vary_L ? length(Ls) : length(ps)
        # data["Mt"] = zeros(nsteps,T) # ⟨magnetization(t)⟩
        scalar_quantities = ["Ft" "binds" "chis" "Ms" "trels"]
        for key in scalar_quantities data[key] = zeros(nsteps) end

        ### compute things requiring monte carlo averages ### 
        function compute(p,q,samps,L,T)
            # define stuff that will be re-used many times by update_state!
            interactions = zeros(L,L)
            for i in 1:L 
                for j in 1:L 
                    interactions[i,j] = 1/circle_distance(i,j,L)^alpha 
                end 
            end
            correction_field = zeros(L)
            correction = zeros(Bool,L)
            charge_correction = zeros(Int,L)

            """ 
            computes monte carlo averages of various quantities for a given value of (p,q,L)
            returns: 
                mode = trel: relaxation time
                mode = erode: final magnetizations 
                mode = Ft: decoding fidelity, time-dependent magnetization
                mode = stats: magnetization, Binder cumulant, susceptibility
            uses synch, eta as global parameters (only p, q, and L may vary)
            """

            mc_keys = ["trels" "Ft" "binds" "chis" "Ms"] # quantities to compute
            mc_data = Dict{String, Any}(key => 0 for key in mc_keys) # store the results of the monte carlo averages
            # mc_data["Mt"] = zeros(T) # time-dependent magnetization

            if mode == "trel"
                maxT = 5000000000
                # logical_mag = 1 # begin the system in the logical state aligned against the bias of the noise 
                @showprogress dt=1 desc="sampling..." for samp in 1:samps 
                    state,charges = charged_init(0,L,charged,alter,"rand") # gives positive magnetization for p < .5
                    t = 1 
                    while t < maxT
                        update_state!(state,charges,charged,p,q,synch,alpha,gamma,interactions,correction_field,correction,charge_correction)
                        if 2*(sum(state)/L - .5) ≤ 0 
                            break 
                        end 
                        t += 1 
                    end
                    if t == maxT 
                        println("maxT reached! (sample = $samp)")
                    end
                    mc_data["trels"] += t / samps 
                end 
            end 

            if mode == "Ft"
                println("T = $T")
                @showprogress dt=1 desc="sampling..." for s in 1:samps 
                    state,charges = charged_init(1,L,"dw",charged) 
                    this_sum = sum(state)
                    
                    for _ in 1:T 
                        # mc_data["Mt"][t] += 2*(this_sum/L-.5) / samps # magnetization 
                        update_state!(state,charges,charged,p,q,synch,alpha,gamma,interactions,correction_field,correction,charge_correction)
                        this_sum = sum(state)
                    end 
                    mc_data["Ft"] += this_sum > L/2 ? 1/samps : 0 # decoding fidelity at end time
                end 
            end

            return mc_data
        end 

        for i in 1:nsteps 

            thisp = p; thisq = q; thisL = L; thissamps = samps[i]
            thisT = Ts[1]
            if vary_L   
                thisL = Ls[i]; thisT = Ts[i]
                println("L = $thisL | samps : $thissamps")
            else
                thisp = ps[i]; thisq = qs[i]
                println("(p,q) = ($thisp,$thisq) | samps : $thissamps")
            end
            this_mc_data = compute(thisp,thisq,thissamps,thisL,thisT)
            for key in keys(this_mc_data)
                data[key][i] = this_mc_data[key]
            end
        end 
    end # done with all calculations 

    # write to file 
    qadj = ""; sadj = ~synch ? "_asynch" : ""; alter_adj = (alter && charged) ? "_alt" : ""; charged_adj = charged ? "charged" : ""; gammaadj = gamma == 0 ? "_gam0" : ""
    if measurement_error_strength > 0 # add flag if measurement noise present 
        qadj = "_q$measurement_error_strength"
    end 
    
    fout = "tmp_data/$(charged_adj)powerdw_$(mode)_alpha$(alpha)_"*(vary_L ? "p$p" : "L$L")*"$sadj$gammaadj$alter_adj$out_adj.jld2"

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