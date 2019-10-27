using Revise
using Random
using SpecialFunctions
using PyCall
using PyPlot
using StatsPlots

function mySliceSampler(pdf_log, x_0, w, m, N)
    # REQUIRES MODULES:
    # Random
    #
    # INPUT:
    # pdf_log --> log-pdf of the target distribution (a function)
    # x_init --> inital vector of x (dimension 1xD), where D dimension of pdf
    # w --> typical slice size, for a window (dimension 1xD)
    # m --> integer limiting the slice size to m*w
    # N --> number of sample points
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


    for ii in 1:N
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
        x_s[ii,:] = x_1    # sample point, after one round of permutations
        pdflog_x_s[ii] = pdflog_x_1 # now contains log_pdf of all new coordinates of new point, since has been updated at each permutation
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

# Jespers HDI:
function hdi(theta_samp,alpha=0.05)
    cred_mass = 1.0-alpha
    ci = zeros(2)
    if length(size(theta_samp))>1   # multidimensional theta
        K,N = size(theta_samp)
        cis = zeros(2,K)
        for k in 1:K    # for each dimension

            ts = theta_samp[k,:]
            sind = sortperm(ts)     # get series of indices that sorts theta values
            sts = ts[sind]  # sorted theta values

            # Shifting a 95% index block from far left to far right and use the
            # samllest interval that can be found during the shift
            N = length(sind)    # number of sample points
            length_ci = Inf     # start with inf HDI interval and shrink
            for i in 1:Int(floor(N*alpha))  # iterate over the 5% lowest values (low to high)
                i2 = Int(floor(N*cred_mass)+i)  # the 5% highest values (low to high), one 95% block away on the other side
                prop_ci = [sts[i],sts[i2]]      # proposed interval
                length_prop_ci = prop_ci[2]-prop_ci[1]
                if length_prop_ci < length_ci
                    ci = prop_ci
                    length_ci = ci[2]-ci[1]
                end
            end
            cis[:,k] = ci

        end
        return cis
    else
        N = length(theta_samp)

        ts = theta_samp
        sind = sortperm(ts)
        sts = ts[sind]

        N = length(sind)
        length_ci = Inf
        for i in 1:Int(floor(N*alpha))
            i2 = Int(floor(N*cred_mass)+i)
            prop_ci = [sts[i],sts[i2]]
            length_prop_ci = prop_ci[2]-prop_ci[1]
            if length_prop_ci < length_ci
                ci = prop_ci
                length_ci = ci[2]-ci[1]
            end
        end
        return ci
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

### Part A
# Chain inputs:
N = 10^6
x_start = 0.5;
w = .5;  # typical window width
m = 100; # maximum number of window widths

#%%
#First sampling
a = 1   # Flat Prior
b = 1
y = [ones(Int64,11); zeros(Int64,3)];   # Given observaitons

log_prior(θ) = log_beta(θ,a,b)
log_lklhd(θ) = log_likelihood_fct(θ,y)
log_posterior(θ) = log_lklhd(θ) + log_prior(θ)

# Need to use "StatsPlots" or "Plots" histogram to add the density curve, PyPlot won't work
chain, pdf = mySliceSampler(log_posterior,x_start,w,m,N);
chain = reshape(chain,:)
#%%
histogram(chain, bins=100, normalize=:pdf, label="posterior")   # Comes from StatsPlots now
density!(chain,linewidth=3,label="posterior")

#%%
# Sample Mean
μ = sum(chain)/N

# Mode
ω = chain[argmax(pdf)]

# Equal Tail Interval:
# Find the 2.5% highest and lowest (equal tails)
# get their index and receive an interval
α = 0.05;
idcesForOrder = sortperm(chain) # Indices of the chain which make it sorted
chain_sorted = chain[idcesForOrder]

leftTailIdx = Int(floor(N*(α/2)))
rightTailIdx = Int(ceil(N*(1- α/2)))

ETI = [chain_sorted[leftTailIdx],chain_sorted[rightTailIdx]]

HDI = hdi(chain)

#%%
p_gthalf = count(i->(i>0.5), chain)/N

#%%

### Part B

# Chain inputs:
N = 10^6
x_start = 0.5;
w = .5;  # typical window width
m = 100; # maximum number of window widths

#%%
#Second dataset
a = 1   # Flat Prior
b = 1
y = [ones(Int64,3); zeros(Int64,7)];   # Given observtions

log_prior(θ) = log_beta(θ,a,b)
log_lklhd(θ) = log_likelihood_fct(θ,y)
log_posterior(θ) = log_lklhd(θ) + log_prior(θ)

# Need to use "StatsPlots" or "Plots" histogram to add the density curve, PyPlot won't work
chain2, pdf2 = mySliceSampler(log_posterior,x_start,w,m,N);
chain2 = reshape(chain2,:);
#%%
histogram(chain2, bins=100, normalize=:pdf, label="posterior")   # Comes from StatsPlots now
density!(chain2,linewidth=3,label="posterior")

#%%  Method 1
gt = (chain .> chain2)*1
sum(gt)/N

#%% Method 1.1, questionable...
#gt = Array{Float64,1}(undef,N)
#for i in 1:length(chain)
    # for each entry in chain 1, compare it to ALL entries in chain2
#    gt[i] = sum((chain[i] .> chain2) * 1)/length(chain2)
#end

#sum(gt)/length(chain)

#%% Method 1.1, questionable...
#gt = (μ .> chain2) * 1
#sum(gt)/N

#%% Method 1.2, questionable...
#μ2 = sum(chain2)/N
#ω2 = chain[argmax(pdf)]

#idcesForOrder = sortperm(chain2) # Indices of the chain which make it sorted
#chain2_sorted = chain[idcesForOrder]

#gt = (chain_sorted .> maximum(chain2_sorted))*1
#sum(gt)

#%% Method 2
dθ = (chain - chain2)
gtZero = (dθ .> 0)*1
sum(gtZero)/N
hdi(dθ)

#%% Method 3
dθ = (chain - chain2)
dθ = reshape(dθ,:)      # to get an Array{Float64,1}
histogram(dθ,bins=100)
diff_hdi = hdi(dθ)
using Plots
Plots.savefig("/home/johhub/Desktop/ABDA/A4/A4-B-hist.pdf")
