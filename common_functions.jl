using Random 
using Statistics
using LinearAlgebra 
using JLD2 
using Alert
using ProgressMeter
using Dates

function circle_distance(i,j,L)
    """
    computes the distance between two sites on a circle of length L
    """
    return min(abs(i-j),L-abs(i-j))
end

function init(p,L,init_type)
    """
    initializes 1d system 
    """
    ind(i) = mod1(i,L)
    state = trues(L)
    if init_type == "rand"
        state = rand(L) .> p
    elseif init_type == "dw" 
        mspins = round(Int,L * p) 
        hL = L ÷ 2
        dx = round(Int,mspins/2)
        state[hL-dx:hL+dx] .= false
    else 
        println("init type $init_type not recognized")
    end 

    return state
end 

function init_2d(p,L,init_type,bconds)
    """
    initializes 2d system
    """
    ind(i) = mod1(i,L)
    state = falses(L,L,2)
    if init_type == "rand"
        if bconds == "periodic"
            for i in 1:L, j in 1:L, o in 1:2
                if rand() < p 
                    state[i,j,o] ⊻= true 
                end 
            end
        else 
            for i in 2:L-1, j in 2:L-1
                if i == L-1 && j < L-1
                    if rand() < p 
                        state[i,j,2] ⊻= true 
                    end 
                end 
                if j == L-1 && i < L-1
                    if rand() < p 
                        state[i,j,1] ⊻= true 
                    end 
                end
                if i < L-1 && j < L-1 
                    for o in 1:2 
                        if rand() < p 
                            state[i,j,o] ⊻= true 
                        end 
                    end
                end 
            end
        end 
    elseif init_type == "dw" 
        qL = L ÷ 4; hL = L ÷ 2
        state[qL:3qL,hL,1] .= false
    else 
        println("init type $init_type not recognized")
    end 

    return state, get_synds(state)
end 

function charged_init(p,L,charged,alter,init_type)
    """
    initializes the state of the system 
    """
    ind(i) = mod1(i,L)
    state = trues(L); charges = zeros(Int,L)
    if init_type == "rand"
        if charged && ~alter  
            for i in shuffle(1:L) 
                if rand() < p 
                    sign = rand([-1,1])
                    if charges[i] != sign && charges[ind(i+1)] != -sign 
                        state[ind(i+1)] ⊻= true 
                        charges[i] += sign 
                        charges[ind(i+1)] -= sign 
                    elseif charges[i] != -sign && charges[ind(i+1)] != sign 
                        state[ind(i+1)] ⊻= true 
                        charges[i] -= sign 
                        charges[ind(i+1)] += sign 
                    end
                end 
            end
        else 
            state = rand(L) .> p
            sign = 1 
            for i in 1:L 
                if state[i] != state[ind(i+1)]
                    charges[i] = sign; sign *= -1 
                end 
            end
            if ~charged 
                charges = abs.(charges)
            end 
        end 
    elseif init_type == "dw" 
        mspins = round(Int,L * p) 
        hL = L ÷ 2
        dx = round(Int,mspins/2)
        state[hL-dx:hL+dx] .= false
        charges[hL-dx-1] = 1 
        charges[hL+dx] = -1 
        if ~charged 
            charges = abs.(charges)
        end 
    else 
        println("init type $init_type not recognized")
    end 

    if charged @assert sum(charges) == 0 "not charge neutral!" end

    return state, charges
end 

function get_synds(state)
    """ 
    calculates syndromes 
    """
    L = size(state)[1]
    ind(i) = mod1(i,L)
    synds = falses(L,L)
    for i in 1:L, j in 1:L
        im1 = ind(i-1); jm1 = ind(j-1)
        synds[i,j] = state[i,j,1] ⊻ state[i,j,2] ⊻ state[im1,j,1] ⊻ state[i,jm1,2]
    end 

    return synds

end 

function detect_logical_error(state)
    """
    input: anyon-free state 
    output: true if both cycles have trivial winding; false otherwise 
    """
    Lx,Ly,_ = size(state)
    xparity = false; yparity = false
    for i in 1:Lx 
        xparity ⊻= state[i,1,2]
    end 
    for j in 1:Ly
        yparity ⊻= state[1,j,1]
    end
    return ~xparity && ~yparity 
end  
