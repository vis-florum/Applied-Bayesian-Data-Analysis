using Revise
using Random
using SpecialFunctions
using PyCall
using PyPlot
using StatsPlots

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
        i = 1 + (ii-1)%N    # gives 1,2,3,...N,1,2,3,...,(N+burnIn)
        x_s[i,:] = x_1    # sample point, after one round of permutations
        pdflog_x_s[i] = pdflog_x_1 # now contains log_pdf of all new coordinates of new point, since has been updated at each permutation
        # note that x_0 has been updated during the permutations

    end

    return x_s, pdflog_x_s
end

#%%
################################################################

function log_beta(θ,a::Any,b::Any)
    #using SpecialFunctions
    θ = θ[1]
    if 0 < θ < 1
        # Numerator
        num = log(θ)*(a-1) + log(1-θ)*(b-1)
        # The denominator (Beta-function) can be integrated or made up of Gamma-functions
        # Since we only treat integer values here, Γ(a) = (a-1)!
        #den = log(factorial(a-1) * factorial(b-1) / factorial(a+b-1))
        den = log(SpecialFunctions.beta(a,b))
        return num - den
    else
        return -Inf
    end
end

################################################################

# Log-Bernoulli PDF:
function log_bernoulli_pdf(θ,y)
    if 0 < θ < 1 && (y==0 || y==1)    # Only admissible values
        return y * log(θ) + (1-y) * log(1-θ)
    else
        return -Inf
    end
end

# Log-Likelihood fct of a (series of) coin tosses with bias:
function log_likelihood_fct(θ, y)
    # Make array s.t. multiple theta values can be input:
    #likelihood = Array{Float64}(undef,length(θ))
    # Do that for every theta value requested:
    #for ii in 1:length(θ)
    #    # Elementwise application of the log_beronoulli_pdf over outcome vector y, then SUM together:
    #    likelihood[ii] = sum(log_bernoulli_pdf.(θ[ii],y))
    #end

    # Un array-ise the output:
    #if length(θ) == 1
    #    likelihood = likelihood[1]
    #end

    θ = θ[1]
    likelihood = sum(log_bernoulli_pdf.(θ,y))

    return likelihood
end

#%%
################################################################

# Chain inputs:
N = 10^6
x_start = 0.5;
w = .5;  # typical window width
m = 100; # maximum number of window widths

# Final posteriors in one array:
posteriors = Array{Float64}(undef,N,3);

#%%
# First Column of Figure 6.4
a = 100
b = 100
y = [ones(1,17) zeros(1,3)]

log_prior(θ) = log_beta(θ,a,b)
log_lklhd(θ) = log_likelihood_fct(θ,y)
log_posterior(θ) = log_lklhd(θ) + log_prior(θ)

# Need to use "StatsPlots" or "Plots" histogram to add the density curve, PyPlot won't work
chain, pdf = mySliceSampler(log_prior,x_start,w,m,N);
histogram(chain, bins=100, normalize=:pdf, label="prior")   # Comes from StatsPlots now
density!(chain,linewidth=3,label="prior")

chain, pdf = mySliceSampler(log_lklhd,x_start,w,m,N)
histogram!(chain, bins=100, normalize=:pdf,label="likelihood")
density!(chain,linewidth=3,label="likelihood")

chain, pdf = mySliceSampler(log_posterior,x_start,w,m,N)
histogram!(chain, bins=100, normalize=:pdf,label="posterior")
density!(chain,linewidth=3,label="posterior")

posteriors[:,1] = chain # save the posterior

Plots.savefig("A4-A-C1_Slice.pdf")

#%%
# Second Column
a = 18.25
b = 6.75
y = [ones(1,17) zeros(1,3)]

log_prior(θ) = log_beta(θ,a,b)
log_lklhd(θ) = log_likelihood_fct(θ,y)
log_posterior(θ) = log_lklhd(θ) + log_prior(θ)

chain, pdf = mySliceSampler(log_prior,x_start,w,m,N);
histogram(chain, bins=100, normalize=:pdf, label="prior")   # Comes from StatsPlots now
density!(chain,linewidth=3,label="prior")

chain, pdf = mySliceSampler(log_lklhd,x_start,w,m,N)
histogram!(chain, bins=100, normalize=:pdf,label="likelihood")
density!(chain,linewidth=3,label="likelihood")

chain, pdf = mySliceSampler(log_posterior,x_start,w,m,N)
histogram!(chain, bins=100, normalize=:pdf,label="posterior")
density!(chain,linewidth=3,label="posterior")

posteriors[:,2] = chain # save the posterior

Plots.savefig("A4-A-C2_Slice.pdf")

#%%
# Third Column
a = 1
b = 1
y = [ones(1,17) zeros(1,3)]

log_prior(θ) = log_beta(θ,a,b)
log_lklhd(θ) = log_likelihood_fct(θ,y)
log_posterior(θ) = log_lklhd(θ) + log_prior(θ)

chain, pdf = mySliceSampler(log_prior,x_start,w,m,N);
histogram(chain, bins=100, normalize=:pdf, label="prior")   # Comes from StatsPlots now
density!(chain,linewidth=3,label="prior")

chain, pdf = mySliceSampler(log_lklhd,x_start,w,m,N)
histogram!(chain, bins=100, normalize=:pdf,label="likelihood")
density!(chain,linewidth=3,label="likelihood")

chain, pdf = mySliceSampler(log_posterior,x_start,w,m,N)
histogram!(chain, bins=100, normalize=:pdf,label="posterior")
density!(chain,linewidth=3,label="posterior")

posteriors[:,3] = chain # save the posterior

Plots.savefig("A4-A-C3_Slice.pdf")

#%%

density(posteriors, label=["a=100, b=100","a=18.25, b=6.75","a=1, b=1"])
Plots.savefig("A4-A-all_Slice.pdf")
