#using PyPlot
#pygui(true)
using Random
using LinearAlgebra
using Plots

# adaptation form Jesper Martinssons implementation of the original slice mySliceSampler:
# https://github.com/jespermartinsson/ABDA.jl/blob/master/src/stats.jl
function mySliceSampler(pdf_log, x_0, w, m, N, burnIn = 1000)
    # REQUIRES MODULES:
    # Random
    #
    # INPUT:
    # pdf_log --> log-pdf of the target distribution (a function)
    # x_init --> inital vector of x (dimension 1xD), where D dimension of pdf
    # w --> typical slice size, for a window (dimension 1xD)
    # m --> integer limiting the slice size to m*w
    # N --> number of sample points
    # burnIn --> number of burn-in samples (optional, results will be discarded)
    #
    # OUTPUT:
    # x_s --> sampled x values (NxD), only keep values after full permutations along dimensions
    # pdflog_x_s --> log-pdf of sampled values

    D = length(x_0)  # the dimension of the distribution (nr of coordinates)
    x_s = Array{Float64}(undef,N,D)   # Samples (NxD), preallocate array, with undefined values
    pdflog_x_s = Array{Float64}(undef,N)   # (unnormalised) log-prob of sampled points (Nx1)
                                           # for reducing evaluations
    pdflog_x_1 = 0.0;

    # to be able to input scalars:
    if typeof(x_0) == Float64
        x_0 = [x_0]
    end


    for ii in 1:(N+burnIn)
        # Update the new x progressively for each dimension:
        L   = 1*x_0  # R and L need to be of same dimension, because we apply the pdf upon those!
        R   = 1*x_0
        x_1 = 1*x_0  # "1*" needed, otherwise x_0 suddendly takes mysterious values
        for d in randperm(D)
            ### 1) Make the Slice
            # random "vertical" position
            if (ii == 1)
                #vertical = rand() * pdf_log(x_0) # in normal space
                z = pdf_log(x_0) + log(rand()) # in log-space
            else
                # saving evaluations, pdflog_x_1 was updated for current permutation in the end
                z = pdflog_x_1 + log(rand())
            end

            ### 2) Make the window, "Stepping-out" alogrithm
            # Randomly place a window, containing x_0:
            L[d] = x_0[d] - rand() * w[d]
            R[d] = L[d] + w[d]

            # Randomly share the max window size among left/right:
            J = floor(m*rand());
            K = (m-1) - J;

            # Extend window to the left, until outside slice or allowance seizes:
            while ((J > 0) && (z < pdf_log(L)))
                #println("Lefting")
                L[d] -= w[d]
                J -= 1
            end

            # Extend window to the right, until outside slice or allowance seizes:
            while ((K > 0) && (z < pdf_log(R)))
                #println("Righting")
                R[d] += w[d]
                K -= 1
            end

            ### 3) Sample from window
            # finding an allowable point:
            while true
                x_1[d] = L[d] + rand() * (R[d] - L[d])

                # this + breaking out of loop reduces the amount of log-pdf evaluations:
                pdflog_x_1 = pdf_log(x_1)

                if (pdflog_x_1 >= z)     # new value found
                    x_0[d] = x_1[d]     # update value for this dimension
                    break
                # Value was not within slice, shrink the interval:
                elseif (x_1[d] < x_0[d])
                    L[d] = x_1[d]
                elseif (x_1[d] > x_0[d])
                    R[d] = x_1[d]
                else
                    throw(ErrorException("Error during shrinking the interval"))
                end
            end

        end # end of permutating through dimensions

        ### 4) Update the chain:
        # Just overwrite the burnIn by using modular arithmetic:
        i = 1 + (ii-1)%N    # gives 1,2,3,...N,1,2,3,...,burnIn
        x_s[i,:] = x_1    # sample point, after one round of permutations
        pdflog_x_s[i] = pdflog_x_1 # now contains log_pdf of all new coordinates of new point, since has been updated at each permutation
        # note that x_0 has been updated during the permutations

    end

    return x_s, pdflog_x_s
end

#%%

using LinearAlgebra

function mylogNormal(x)
    x = x[1]
    m = 1;
    s = 1;
    y = log(1/sqrt(2*pi*s^2)) + (-(x-m)^2 / (2*s^2))
    return y
end

function mylogNormal_bivariate(x)
    # Requires LinearAlgebra
    X = [0.0; 0.0]
    X[1] = x[1]
    X[2] = x[2]

    M = [1;
         1]
    S = [1 0;
         0 1]

    y = log(1 / sqrt(2*pi)) +
        log(1 / sqrt(det(S))) +
        (-1/2) * (X-M)' * inv(S) * (X-M)

    return y
end


function myCrazyDistribution(x)
    if x < 1
        y = 0
    elseif 1 <= x <= 2
        y = 1
    elseif 2 < x < 3
        y = 0
    elseif 3 <= x <= 4
        y = 1
    else
        y = 0
    end
    return y
end


function log_square(x)
    x=x[1]
    if 0 < x < 4
        y = 2*log(x)
    else
        y = -Inf
    end
    return y
end

function log_uniform(x)
    x = x[1]
    if 0 < x < 1
        y = 0
    else
        y = -Inf
    end
end


#%%
N = 1000000
#x_start = [0.5];
x_start = [0.5 0.5];
w = [.1 0.1]; # typical window width
m = 100;

using Distributions
bivarLogNormal(X) = log(pdf(MvNormal([1.0;1.0],[1. 0.; 0. 1.]), [X[1];X[2]]))

bivarLogNormal([1 1])

#%%
#x_s, log_vals = mySliceSampler(log_uniform, x_start, w, m, N)
#x_s, log_vals = mySliceSampler(log_square, x_start, w, m, N)
#x_s, log_vals = mySliceSampler(mylogNormal, x_start, w, m, N)
#x_s, log_vals = mySliceSampler(bivarLogNormal, x_start, w, m, N)
x_s, log_vals = mySliceSampler(mylogNormal_bivariate, x_start, w, m, N)

mean = sum(x_s,dims=1)/N

#%%
#figure()
#hist(x_s[:,2],bins=100)
histogram2d(x_s[:,1],x_s[:,2])#bins=100)

#%%

x = -2:.1:4
y = -2:.1:4
f(x,y) = exp(mylogNormal_bivariate([x;y]))
pyplot()
plot(x,y,f,st=:surface)


#%%
