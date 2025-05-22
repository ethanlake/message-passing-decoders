using Random 
using Statistics
using LinearAlgebra 
using JLD2 
using Alert
using Plots
using ProgressMeter
using Dates

"""
uses random noise to create a distribution of anyons in 1 or 2 dimensions. keeps track of which anyons are paired with which, and calculates statistics of the distribution of pair distances 
"""

function closest_neighbor(A,i)  
    """
    returns the dist to the closest nonzero element of A to position i 
    """
    L = length(A)
    ind(i) = mod1(i,L)
    hL = L÷2
    dx = 1 
    while dx < hL 
        if A[ind(i+dx)] != 0 || A[ind(i-dx)] != 0
            return dx 
        end
        dx += 1
    end
end

function closest_unpaired_neighbor(A,i)  
    """
    returns the dist to the closest nonzero element of A to position i which is not paired with the element at position i 
    """
    L = length(A)
    ind(i) = mod1(i,L)
    hL = L÷2
    dx = 1 
    avoid = [0, i]
    while dx < hL 
        if ~(A[ind(i+dx)] ∈ avoid)  || ~(A[ind(i-dx)] ∈ avoid)
            return dx 
        end
        dx += 1
    end
    return dx # return half the system size if only have a single pair in the system 
end

function closest_neighbor(A, i, j)
    """
    returns 1-norm dist to the closest nonzero element of A to position i,j
    """
    L = size(A)[1]; hL = L÷2 # will only need to search this far in practice
    dist = 1 
    ind(i) = mod1(i,L)
    while dist < hL
        for s in 0:dist
            println(s)
            if A[ind(i+s),ind(j+dist-s),1] != 0 || A[ind(i-s),ind(j+dist-s),1] != 0 || A[ind(i+s),ind(j-dist+s),1] != 0 || A[ind(i-s),ind(j-dist+s),:] != 0 
                return dist  
            end 
        end 
        dist += 1 
    end
    return dist 
end

function closest_unpaired_neighbor(A, i, j)
    """
    returns 1-norm dist to the closest nonzero element of A to position i,j
    """
    L = size(A)[1]; hL = L÷2 # will only need to search this far in practice
    dist = 1 
    ind(i) = mod1(i,L)
    avoid = [[0,0], [i,j]]
    while dist < hL
        for s in 0:dist
            if ~(A[ind(i+s),ind(j+dist-s),:] ∈ avoid) || ~(A[ind(i-s),ind(j+dist-s),:] ∈ avoid) || ~(A[ind(i+s),ind(j-dist+s),:] ∈ avoid) || ~(A[ind(i-s),ind(j-dist+s),:] ∈ avoid)
                return dist  
            end 
        end 
        dist += 1 
    end
    if dist == hL 
        return hL
    end 
end

function circle_dist(i,j,L)
    """
    returns the dist between two points on a circle with L points
    """
    return min(abs(i-j),L-abs(i-j)) 
end

function torus_dist(r1,r2,L)
    """
    returns the 1-norm dist between two points on a torus with side length L
    """ 
    bigx = max(r1[1],r2[1]); smallx = min(r1[1],r2[1])
    bigy = max(r1[2],r2[2]); smally = min(r1[2],r2[2])

    xdist = min(bigx - smallx,L + smallx - bigx)
    ydist = min(bigy - smally,L + smally - bigy)

    return xdist + ydist
end

function generate_state(p,L,d)
    """
    generates a state on a d = {1,2} dimensional lattice with particles at random positions. noise is iid with strength p. 
    """
    ind(i) = mod1(i,L)
    if d == 1 
        positions = zeros(Int,L)
        for i in 1:L 
            if rand() < p 
                ip1 = ind(i+1)
                if positions[i] == 0 
                    if positions[ip1] == 0 # particle pair created 
                        positions[i] = ip1 
                        positions[ip1] = i
                    else # particle moved from site i+1 to site i 
                        positions[i] = positions[ip1]
                        positions[positions[ip1]] = i
                        positions[ip1] = 0
                    end
                else 
                    if positions[ip1] == 0 # particle moved from site i to site i+1 
                        positions[ip1] = positions[i]
                        positions[positions[i]] = ip1
                        positions[i] = 0
                    else # particle pair destroyed 
                        positions[positions[i]] = positions[ip1]
                        positions[positions[ip1]] = positions[i]
                        positions[i] = 0
                        positions[ip1] = 0
                    end
                end
            end 
        end 
        return positions

    elseif d == 2 # two dimensions

        positions = zeros(Int,L,L,2) # positions[i,j,k] = kth coordinate of particle paired to the particle at i,j if nonzero

        for i in 1:L, j in 1:L, o in 1:2 
            if rand() < p 
                ip1 = ind(i+1); jp1 = ind(j+1)
                if o == 1 # add link along x direction 
                    if positions[i,j,1] == 0 # no particle already present 
                        if positions[ip1,j,1] == 0 # particle pair created 
                            positions[i,j,1] = ip1 
                            positions[i,j,2] = j
                            positions[ip1,j,1] = i
                            positions[ip1,j,2] = j
                        else # particle moved from site i+1,j to site i,j
                            positions[i,j,:] .= positions[ip1,j,:]
                            positions[positions[ip1,j,1],positions[ip1,j,2],1] = i
                            positions[positions[ip1,j,1],positions[ip1,j,2],2] = j
                            positions[ip1,j,:] .= 0
                        end
                    else # particle already present 
                        if positions[ip1,j,1] == 0 # particle moved from site i to site i+1 
                            positions[ip1,j,:] .= positions[i,j,:]
                            positions[positions[i,j,1],positions[i,j,2],1] = ip1
                            positions[positions[i,j,1],positions[i,j,2],2] = j
                            positions[i,j,:] .= 0
                        else # particle pair destroyed 
                            positions[positions[i,j,1],positions[i,j,2],1] = positions[ip1,j,1]
                            positions[positions[i,j,1],positions[i,j,2],2] = positions[ip1,j,2]
                            positions[positions[ip1,j,1],positions[ip1,j,2],1] = positions[i,j,1]
                            positions[positions[ip1,j,1],positions[ip1,j,2],2] = positions[i,j,2]
                            positions[i,j,:] .= 0
                            positions[ip1,j,:] .= 0
                        end
                    end
                else # add link along y direction
                    if positions[i,j,1] == 0 # no particle already present
                        if positions[i,jp1,1] == 0 # particle pair created 
                            positions[i,j,1] = i 
                            positions[i,j,2] = jp1
                            positions[i,jp1,1] = i
                            positions[i,jp1,2] = j
                        else # particle moved from site i,j+1 to site i,j
                            positions[i,j,:] .= positions[i,jp1,:]
                            positions[positions[i,jp1,1],positions[i,jp1,2],1] = i
                            positions[positions[i,jp1,1],positions[i,jp1,2],2] = j
                            positions[i,jp1,:] .= 0
                        end
                    else # particle already present
                        if positions[i,jp1,1] == 0 # particle moved from site j to site j+1 
                            positions[i,jp1,:] .= positions[i,j,:]
                            positions[positions[i,j,1],positions[i,j,2],1] = i
                            positions[positions[i,j,1],positions[i,j,2],2] = jp1
                            positions[i,j,:] .= 0
                        else # particle pair destroyed 
                            positions[positions[i,j,1],positions[i,j,2],1] = positions[i,jp1,1]
                            positions[positions[i,j,1],positions[i,j,2],2] = positions[i,jp1,2]
                            positions[positions[i,jp1,1],positions[i,jp1,2],1] = positions[i,j,1]
                            positions[positions[i,jp1,1],positions[i,jp1,2],2] = positions[i,j,2]
                            positions[i,j,:] .= 0
                            positions[i,jp1,:] .= 0
                        end
                    end
                end 
            end 
        end
        return positions
    else 
        println("higher dimensions not yet implemented")
    end 
end 

function get_stats(positions)
    """
    computes statistics of a particle configuration: average and largest dist between particles in a pair, and density of anyons
    """
    d = length(size(positions)) == 3 ? 2 : 1
    L = size(positions)[1]

    mean_pair_dist = 0  
    max_pair_dist = 0
    mean_neighbor_dist = 0 
    density = 0 
    N = 0 
    nearest_is_pair = 0 

    if d == 1 
        for i in 1:L 
            if positions[i] != 0 
                pair_dist = circle_dist(i,positions[i],L)
                mean_pair_dist += pair_dist 
                unpaired_dist = closest_unpaired_neighbor(positions,i)
                mean_neighbor_dist += unpaired_dist

                if unpaired_dist > pair_dist nearest_is_pair += 1 end 

                if pair_dist > max_pair_dist 
                    max_pair_dist = pair_dist 
                end
                N += 1
            end
        end
        mean_pair_dist /= N; mean_neighbor_dist /= N
        density = N / L
        nearest_is_pair /= N

    elseif d == 2 
        for i in 1:L, j in 1:L 
            if positions[i,j,1] != 0 
                pair_dist = torus_dist((i,j),(positions[i,j,1],positions[i,j,2]),L)
                mean_pair_dist += pair_dist 
                unpaired_dist = closest_unpaired_neighbor(positions,i,j)
                mean_neighbor_dist += unpaired_dist
                if unpaired_dist > pair_dist nearest_is_pair += 1 end

                if pair_dist > max_pair_dist 
                    max_pair_dist = pair_dist 
                end
                N += 1
            end
        end
        mean_pair_dist /= N; mean_neighbor_dist /= N
        density = N / (L^2)
        nearest_is_pair /= N
    end

    return mean_pair_dist, max_pair_dist, mean_neighbor_dist, nearest_is_pair, density 
end 

function main() 

    d = 2
    p = 0.5
    L = d == 1 ? 500 : 50 
    nps = 10 # number of noise steps 
    out_adj = "" # suffix of output file name 

    vary_L = ~true 

    Lmin = 25; Lmax = 2500
    if d == 2 
        Lmin = 8; Lmax = 64
    end

    Ls = vary_L ? [round(Int,2^el) for el in LinRange(log2(Lmin),log2(Lmax),nps)] : [L for _ in 1:nps]

    pmin = .01; pmax = .25
    ps = vary_L ? [p for _ in 1:nps] : [pe for pe in range(pmin,pmax,length=nps)]
    samps = [round(Int,70000 * (d == 1 ? (el / 100) : (el / 25)^2)) for el in Ls]

    means,maxs,densities,mean_nns,nearest_is_pairs = [zeros(nps) for _ in 1:5] # various statistical quantities we want to measure 
    
    for pind in 1:nps 
        thisL = Ls[pind]; thisp = ps[pind]; thissamps = samps[pind]
        println("p = $thisp; L = $thisL")
        @showprogress dt = 1 for _ in 1:thissamps 
            state = generate_state(thisp,thisL,d)
            mean_pair_dist, max_pair_dist, mean_neighbor_dist, nearest_is_pair, density = get_stats(state)
            means[pind] += mean_pair_dist
            maxs[pind] += max_pair_dist
            mean_nns[pind] += mean_neighbor_dist
            nearest_is_pairs[pind] += nearest_is_pair
            densities[pind] += density
        end
        means[pind] /= thissamps; maxs[pind] /= thissamps; densities[pind] /= thissamps; mean_nns[pind] /= thissamps; nearest_is_pairs[pind] /= thissamps
    end 

    params = Dict{String,Any}("d"=>d,"p"=>p,"L"=>L,"nps"=>nps,"vary_L"=>vary_L,"Ls"=>Ls,"ps"=>ps,"samps"=>samps)
    data = Dict{String,Any}("means"=>means,"maxs"=>maxs,"densities"=>densities,"mean_nns"=>mean_nns,"nearest_is_pairs"=>nearest_is_pairs)

    fout = "data/"*(d==2 ? "2d" : "1d")*"pair_dists_"*(vary_L ? "p$p" : "L$L")*out_adj*".jld2"
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