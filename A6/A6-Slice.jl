#using Revise
using Distributed
@everywhere using Random
@everywhere using Statistics
using Plots
using StatsPlots
using SharedArrays


################################################################################
############################# PREPARATION ######################################

#%% Input Data #################################################################
# Individual observations in sequence:
@everywhere y = [607, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395, 346, 361, 310,
300, 382, 294, 315, 323, 421, 339, 398, 328, 335, 291, 329, 310, 294, 321, 286,
349, 279, 268, 293, 310, 259, 241, 243, 272, 247, 275, 220, 245, 268, 357, 273,
301, 322, 276, 401, 368, 149, 507, 411, 362, 358, 355, 362, 324, 332, 268, 259,
274, 248, 254, 242, 286, 276, 237, 259, 251, 239, 247, 260, 237, 206, 242, 361,
267, 245, 331, 357, 284, 263, 244, 317, 225, 254, 253, 251, 314, 239, 248, 250,
200, 256, 233, 427, 391, 331, 395, 337, 392, 352, 381, 330, 368, 381, 316, 335,
316, 302, 375, 361, 330, 351, 186, 221, 278, 244, 218, 126, 269, 238, 194, 384,
154, 555, 387, 317, 365, 357, 390, 320, 316, 297, 354, 266, 279, 327, 285, 258,
267, 226, 237, 264, 510, 490, 458, 425, 522, 927, 555, 550, 516, 548, 560, 545,
633, 496, 498, 223, 222, 309, 244, 207, 258, 255, 281, 258, 226, 257, 263, 266,
238, 249, 340, 247, 216, 241, 239, 226, 273, 235, 251, 290, 473, 416, 451, 475,
406, 349, 401, 334, 446, 401, 252, 266, 210, 228, 250, 265, 236, 289, 244, 327,
274, 223, 327, 307, 338, 345, 381, 369, 445, 296, 303, 326, 321, 309, 307, 319,
288, 299, 284, 278, 310, 282, 275, 372, 295, 306, 303, 285, 316, 294, 284, 324,
264, 278, 369, 254, 306, 237, 439, 287, 285, 261, 299, 311, 265, 292, 282, 271,
268, 270, 259, 269, 249, 261, 425, 291, 291, 441, 222, 347, 244, 232, 272, 264,
190, 219, 317, 232, 256, 185, 210, 213, 202, 226, 250, 238, 252, 233, 221, 220,
287, 267, 264, 273, 304, 294, 236, 200, 219, 276, 287, 365, 438, 420, 396, 359,
405, 397, 383, 360, 387, 429, 358, 459, 371, 368, 452, 358, 371];

# Indicator of the individual behind each observation
@everywhere ind = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5,
5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8,
8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14,
14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16,
16, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 20, 20, 20,
20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21,
21, 21, 21, 21, 21, 21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22, 22,
22, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24,
24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25,
25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 28, 28, 28,
28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30,
30, 30, 30, 31, 31, 31, 31, 31, 32, 32, 32, 32, 32, 33, 34, 34, 34, 34, 34,
34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34];

# Indicator for each individual j whether he/she is a child or not:
@everywhere child_j = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1,
0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];

# Indicator for each observations i whether it comes from a child or not:
@everywhere child_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1,
1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];


#####
#%% Functions ##################################################################
# @everywhere to make them known to all workers
@everywhere function mySliceSampler(pdf_log, x_0, w, m, N = 2000, burnIn = 1000)
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
    #if typeof(x_0) == Float64
    #    x_0 = [x_0]
    #end


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


## Jespers HDI:
@everywhere function hdi(theta_samp,alpha=0.05)
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


# A recursive application of the HDI algorithm to find the mode:
@everywhere function find_mode(sample_pts;runs=2)
    # re-apply the hdi algorithm, until there is only two points or so left
    # so the ci interval becomes smaller and smaller
    # in the end, return an average over the interval
    ci = zeros(2)
    i_left = 0
    i_right = 0

    # Make 50% credibility intervals
    alpha = 0.5;
    cred_mass = 1.0 - alpha;

    # sort the pts:
    sind = sortperm(sample_pts)
    sts = sample_pts[sind]
    N = length(sind)

    if N > runs
        # HDI algorithm, but return indices too:
        length_ci = Inf
        for i in 1:Int(floor(N*alpha))
            i2 = Int(floor(N*cred_mass)+i)
            prop_ci = [sts[i], sts[i2]]
            length_prop_ci = prop_ci[2] - prop_ci[1]
            if length_prop_ci < length_ci
                i_left  = i
                i_right = i2
                ci = prop_ci
                length_ci = ci[2] - ci[1]
            end
        end

        # Recursion:
        ci = find_mode(sts[i_left:i_right],runs=runs)

    else
        ci = [sts[1],sts[end]]
    end

    return ci
end


# A function to make the plots with HDIs
function makeDistributionPlot(X, color="blue"; ann=true, offset=0.0, scale=1.0)
    #%% Mode, Mean, HDIs
    ω = mean(find_mode(X))
    μ_bar = mean(X);
    left,right = hdi(X);

    h = fit(Histogram, X, nbins=100)
    w = h.weights      # weights of each bar
    e = h.edges[1]     # edges of the bars (must be +1 more than bars)

    riemannSum = sum(w.*diff(e))
    w = w ./ riemannSum

    e2 = similar([e[:]; e[:]])  # doubling of the edge vector
    e2[2:2:end] = e             # even indices
    e2[1:2:end] = e             # odd indices
    e2 = e2[2:(end-1)]          # pop the ends

    w2 = similar([w[:]; w[:]])  # doubling of the weight vector
    w2[2:2:end] = w             # even indices
    w2[1:2:end] = w             # odd indices
    w2 = w2 .* scale .+ offset

    w_left = w[sum(e .< left)] * scale + offset  # comparison leaves at max left edge of bar that contains the limit
    w_right = w[sum(e .< right)] * scale + offset
    w_ω = w[sum(e .< ω)] * scale + offset
    w_μ = w[sum(e .< μ_bar)] * scale + offset

    #%% Plotting
    plot(e2, w2,
         linealpha=0.1, linecolor=color,
         fillrange=offset, fillalpha=0.4, fillcolor=color,
         legend=false)

    # or do it by histogram, then we also have vertical lines
    # specify statsplot?
    #histogram(X, bins=100, normalize=:pdf,
    #          color=color, alpha=0.4, linealpha=0.1, legend=false)   # Comes from StatsPlots now

    if ann
    plot!([(left,offset),(left,w_left)],
        linewidth=1.5, color=color, linestyle=:dash,
        annotations = (left, w_left, text("$(round(left,sigdigits=4))",color,rotation=90,:center,:left)))
    plot!([(right,offset),(right,w_right)],
        linewidth=1.5, color=color, linestyle=:dash,
        annotations = (right, w_right, text("$(round(right,sigdigits=4))",color,rotation=90,:center,:left)))
    plot!([(μ_bar,offset),(μ_bar,w_μ)],
        linewidth=1.5, color=color, linestyle=:dash,
        annotations = (μ_bar, 0.5*w_μ, text("mean $(round(μ_bar,sigdigits=4))",color,rotation=90,:bottom,:right)))
    plot!([(ω,offset),(ω,w_ω)],
        linewidth=1.5, color=color, linestyle=:dashdot,
        annotations = (ω, w_ω, text("mode $(round(ω,sigdigits=4))",color,rotation=90,:bottom,:right)))
    else
    plot!([(left,offset),(left,w_left)],
        linewidth=1.5, color=color, linestyle=:dash)
    plot!([(right,offset),(right,w_right)],
        linewidth=1.5, color=color, linestyle=:dash)
    plot!([(μ_bar,offset),(μ_bar,w_μ)],
        linewidth=1.5, color=color, linestyle=:dash)
    plot!([(ω,offset),(ω,w_ω)],
        linewidth=1.5, color=color, linestyle=:dashdot)
    end
end


function makeDistributionPlot!(X, color="blue"; ann=true, offset=0.0, scale=1.0)
    #%% Mode, Mean, HDIs
    ω = mean(find_mode(X))
    μ_bar = mean(X);
    left,right = hdi(X);

    h = fit(Histogram, X, nbins=100)
    w = h.weights      # weights of each bar
    e = h.edges[1]     # edges of the bars (must be +1 more than bars)

    riemannSum = sum(w.*diff(e))
    w = w ./ riemannSum

    e2 = similar([e[:]; e[:]])  # doubling of the edge vector
    e2[2:2:end] = e             # even indices
    e2[1:2:end] = e             # odd indices
    e2 = e2[2:(end-1)]          # pop the ends

    w2 = similar([w[:]; w[:]])  # doubling of the weight vector
    w2[2:2:end] = w             # even indices
    w2[1:2:end] = w             # odd indices
    w2 = w2 .* scale .+ offset

    w_left = w[sum(e .< left)] * scale + offset  # comparison leaves at max left edge of bar that contains the limit
    w_right = w[sum(e .< right)] * scale + offset
    w_ω = w[sum(e .< ω)] * scale + offset
    w_μ = w[sum(e .< μ_bar)] * scale + offset

    #%% Plotting
    plot!(e2, w2,
         linealpha=0.1, linecolor=color,
         fillrange=offset, fillalpha=0.4, fillcolor=color,
         legend=false)

    # or do it by histogram, then we also have vertical lines
    # specify statsplot?
    #histogram!(X, bins=100, normalize=:pdf,
    #          color=color, alpha=0.4, linealpha=0.1, legend=false)   # Comes from StatsPlots now

    if ann
    plot!([(left,offset),(left,w_left)],
        linewidth=1.5, color=color, linestyle=:dash,
        annotations = (left, w_left, text("$(round(left,sigdigits=4))",color,rotation=90,:center,:left)))
    plot!([(right,offset),(right,w_right)],
        linewidth=1.5, color=color, linestyle=:dash,
        annotations = (right, w_right, text("$(round(right,sigdigits=4))",color,rotation=90,:center,:left)))
    plot!([(μ_bar,offset),(μ_bar,w_μ)],
        linewidth=1.5, color=color, linestyle=:dash,
        annotations = (μ_bar, 0.5*w_μ, text("mean $(round(μ_bar,sigdigits=4))",color,rotation=90,:bottom,:right)))
    plot!([(ω,offset),(ω,w_ω)],
        linewidth=1.5, color=color, linestyle=:dashdot,
        annotations = (ω, w_ω, text("mode $(round(ω,sigdigits=4))",color,rotation=90,:bottom,:right)))
    else
    plot!([(left,offset),(left,w_left)],
        linewidth=1.5, color=color, linestyle=:dash)
    plot!([(right,offset),(right,w_right)],
        linewidth=1.5, color=color, linestyle=:dash)
    plot!([(μ_bar,offset),(μ_bar,w_μ)],
        linewidth=1.5, color=color, linestyle=:dash)
    plot!([(ω,offset),(ω,w_ω)],
        linewidth=1.5, color=color, linestyle=:dashdot)
    end
end

# Taken from Jesper:
# https://en.wikipedia.org/wiki/Correlation_coefficient
function ACF(x, k=1)
  x_meanshift = x .- mean(x)
  zxk = x_meanshift[k+1:end]
  zyk = x_meanshift[1:end-k]
  return sum(zxk.*zyk)/sqrt(sum(zxk.^2)*sum(zyk.^2))
end

# Taken from Jesper:
function acovlim(x;lim=0.05)
  k = 0
  rhos = []
  rho = 1
  while rho>lim
    rho = ACF(x,k)
    push!(rhos,rho)
    k += 1
  end
  return rhos
end

# Taken from Jesper:
# ess -- effective sample size (Kruschke 2014, page 184)
function ess(x)
    if typeof(x)==Vector{Float64}
        n = length(x)
        acf = acovlim(x)
        return n/(1+2*sum(acf[2:end]))
    else
        m,n = size(x)
        list = zeros(m)
        for i in 1:m
            acf = acovlim(x[i,:])
            list[i] = n/(1+2*sum(acf[2:end]))
        end
        return list
    end
end


################################################################################

#####
#%% Priors, Likelihoods, Posteriors ############################################
#####
@everywhere function log_prior_flat(X::Float64)
    return 0.0  # or any other constant value, for the posterior it doesn't matter
end


@everywhere function log_prior_flat_positiveOnly(X::Float64)
    if X > 0
        return 0.0  # or any other constant value, for the posterior it doesn't matter
    else
        return -Inf
    end
    # Alternative of Jesper:
    #return log(tau > 0)
end


# Here we have to
# accomodate, that the means no longer are a constant, but an array, i.e. a
# a mean for each theta
@everywhere function log_prior_theta(θ::Array{Float64},μ::Array{Float64},σ::Float64)
    # Log-normal distribution
    # Logarithmise the normal pdfs in series:
    # We can leave out the 1/sqrt(2π) since it does not change sampling from the posterior
    return sum(-log(σ) .- 0.5 .* ((θ .- μ) ./ σ).^2)
end

@everywhere function log_lklhd_fct(jointParams::Array{Float64},zlogy::Array{Float64,1},ind::Array{Int64,1})
    # return the joint lklhd of the given input data
    J = length(jointParams) - 4
    θ = jointParams[1:J]
    μ_0, σ, τ, ϕ  = jointParams[(J+1):end]

    jointLklhd = 0.0

    if (σ > 0)
        # Logarithmise the normal pdfs in series:
        # summing up:
        # note, that we go through each observation, vectorising is hard:
        for i in 1:length(zlogy)
            jointLklhd += -log(σ) - 0.5 * ((zlogy[i] - θ[ind[i]]) / σ)^2
        end
    else
        jointLklhd = -Inf
    end

    return jointLklhd
end


@everywhere function log_priors(jointParams::Array{Float64},K_j::Array{Int64,1})
    # The joint log-pdf of all the priors in the model
    J = length(jointParams) - 4
    θ = jointParams[1:J]
    μ_0, σ, τ, ϕ  = jointParams[(J+1):end]

    μ = μ_0 .+ ϕ .* K_j.*1.0    # current sampled μ for each individual

    if (σ > 0) && (τ > 0)
        return log_prior_theta(θ,μ,τ) +
               log_prior_flat(μ_0) +
               log_prior_flat(ϕ) +
               log_prior_flat_positiveOnly(σ) +
               log_prior_flat_positiveOnly(τ)
    else
        return -Inf
    end
end

#%% Logarithmise and Transform Data ############################################
J = maximum(ind)    # number of individuals
I = length(y)       # number of observations

# Go to log-space:
@everywhere logy = log.(y)
@everywhere logMean = mean(logy)
@everywhere logStd = std(logy)

# Mean-centre and scale each element:
@everywhere zlogy = (logy .- logMean) ./ logStd


#%% Concretise the functions by including the data #############################
@everywhere log_lklhd(jointParams::Array{Float64}) =
            log_lklhd_fct(jointParams::Array{Float64}, zlogy, ind)
@everywhere log_posterior(jointParams::Array{Float64}) =
            log_lklhd(jointParams::Array{Float64}) +
            log_priors(jointParams::Array{Float64},child_j)



################################################################################
############################# MCMC #############################################
#%% Chain Setup
noOfChains = 4
N = 5*10^6                # total number of samples
N_chain = convert(Int64, N / noOfChains)   # samples per chain
burnIn = 10^4
w = ones(Float64,J+4) * 0.1         # typical window size
m = 100                             # multiplier for maximum window size

#%% MCMC run, make several chains
chain = SharedArray{Float64}(N_chain,J+4,noOfChains);
chainPdf = SharedArray{Float64}(N_chain,noOfChains);

# (using threads for parallel computing (other stuff didn't work)):
# Threads.@threads led to problems too
# now using @distributed, prefix with @sync to wait for completion
@sync @distributed for k in 1:noOfChains
    x_start = rand(Float64,J+4) * 0.5   # initial vector
    chain[:,:,k], chainPdf[:,k] = mySliceSampler(log_posterior,x_start,w,m,N_chain,burnIn);
end
# My machine: 4 cores, i7-7500U 2,7 GHz x 2 each, on Linux
# N=10^5 takes around 25 seconds.
# N=10^6 takes around 120 seconds.
# N=5*10^6 takes around 590 seconds.


################################################################################
############################# RESULTS ##########################################

##############
# Acess results and Transform back (undo mean-centering and scaling and go to non-log space)
# See derivation in PDF file instead

# Access the parameters and fuse the chains:
θ = Array{Float64,2}(undef,N,J)
for j in 1:J
    θ[:,j] = chain[:,j,:][:]
end
μ_0 = chain[:,J+1,:][:]
σ   = chain[:,J+2,:][:]
τ   = chain[:,J+3,:][:]
ϕ   = chain[:,J+4,:][:]

#####
# Autocorrelation and ESS
ess_θ = Array{Float64,1}(undef,J)
for j in 1:J
    ess_θ[j] = ess(θ[:,j])
    println("ESS of θ_$j = ",ess_θ[j])
end
println("Minimal ESS is of θ_$(argmin(ess_θ)) = ",ess(θ[:,argmin(ess_θ)]))
println("Maximal ESS is of θ_$(argmax(ess_θ)) = ",ess(θ[:,argmax(ess_θ)]))

println("ESS of μ_0 = ",ess(μ_0))
println("ESS of σ = ",ess(σ))
println("ESS of τ = ",ess(τ))
println("ESS of ϕ = ",ess(ϕ))


#####
# Un-scale and un-mean-centre:
θ_unscaled = θ .* logStd .+ logMean
# adults:
μ_0_unscaled = (μ_0 .+ 0) .* logStd .+ logMean
# kids only:
μ_ϕ_unscaled = (μ_0 .+ ϕ) .* logStd .+ logMean
ϕ_unscaled = ϕ .* logStd

σ_unscaled = σ .* logStd
τ_unscaled = τ .* logStd

#####
# Get into non-log space:
θ_unscaled_unLog = exp.(θ_unscaled .+ 0.5 .* repeat(σ_unscaled,1,J).^2);
μ_0_trans_unLog = exp.(μ_0_unscaled .+ 0.5 .* σ_unscaled.^2 .+ 0.5 .* τ_unscaled.^2);    # adults
μ_ϕ_trans_unLog = exp.(μ_0_unscaled .+ ϕ_unscaled .+ 0.5 .* σ_unscaled.^2 .+ 0.5 .* τ_unscaled.^2);  # kids


################################################################################
############################# TASKS ############################################

######################## Task 1 ####################################
# Effect of being a kid:
makeDistributionPlot(ϕ, "orange")
Plots.savefig(projDir*"/figs/phi_Slice.pdf")
makeDistributionPlot(ϕ_unscaled, "orange")
Plots.savefig(projDir*"/figs/phi_unscaled_Slice.pdf")


######################## Task 2 #############
makeDistributionPlot(τ,"blue")
Plots.savefig(projDir*"/figs/tau_Slice.pdf")

makeDistributionPlot(τ_unscaled, "blue")
Plots.savefig(projDir*"/figs/tau_unscaled_Slice.pdf")


######################## Task 3, Priors of expected log reaction #############
# prior for theta was
# theta[j] ~ normal(mu + phi*ISKID[j],tau)

prior_adult = mean(μ_0_unscaled) .+ mean(τ_unscaled) .* randn(N)
prior_kid = mean(μ_0_unscaled) .+ mean(ϕ_unscaled) .+ mean(τ_unscaled) .* randn(N)
makeDistributionPlot(prior_adult,"black")
makeDistributionPlot!(prior_kid,"red")
Plots.savefig(projDir*"/figs/priors_Slice.pdf")

# Logarithmic sample means
θ_mean_unscaled = zeros(J)
for j in 1:J
    θ_mean_unscaled[j] = mean(logy[ind .== j])
end

#%% Plot into one figure
# Get ordered indices:
sIdx = sortperm(θ_mean_unscaled)

plt = makeDistributionPlot(prior_adult,"black",ann=false)
makeDistributionPlot!(prior_kid,"red",ann=false)
top = 1.6;
for i in 1:J
    j = sIdx[i]
    if child_j[j] == 1
        color = "red"
    else
        color = "black"
    end
    makeDistributionPlot!(θ_unscaled[:,j],color,ann=false,offset=i*top/J,scale=1/100)

end
plt

Plots.savefig(projDir*"/figs/priors_postOverlay_Slice.pdf")


######################## Task 4, posterior prediction #############
#### a) knowing that it is a child

# Posterior predictive sampling:
idx = Int.(ceil.(rand(N).*N))   # choose random indices which will pick from the simulated posterior
# Sample zlogy according to the model:
zlogy_sim_adult = μ_0[idx] .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
zlogy_sim_kid = μ_0[idx] .+ ϕ[idx] .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
# Transform to non-standardised and non-log:
y_sim_adult = exp.(zlogy_sim_adult .* logStd .+ logMean)
y_sim_kid = exp.(zlogy_sim_kid .* logStd .+ logMean)

makeDistributionPlot(y_sim_adult,"black",ann=true)
makeDistributionPlot!(y_sim_kid,"red",ann=true)
Plots.savefig(projDir*"/figs/PP_known_Slice.pdf")

##############
#### b) not knowing that it is a child
# assuming a prior of equally likely to be adult or child (beta(1,1))
# then add a likelihood corresponding to a bernoulli processs, with number of
# heads equvivalent to number of kids among all individuals
# the posterior probability of having the amount of childs among
# our total individuals is given by a bernoulli process with a beta
# with parameters a=b=1, N = # of individuals, z = # of kids

nr_kids = sum(child_j)  # number of heads
a = nr_kids + 1         # head (kid) count here
b = J - nr_kids + 1     # non-heads (adults)
postBeingKid = rand(Beta(a,b),N)
areKids = postBeingKid .>= rand(N)  # only those who make it over the threshold

# Posterior predictive sampling
idx = Int.(ceil.(rand(N).*N))   # choose a random index which will pick from the simulated posterior
zlogy_sim_unknown = μ_0[idx] .+ ϕ[idx] .* areKids .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
y_sim_unknown = exp.(zlogy_sim_unknown .* logStd .+ logMean)
# Compare in plots:
makeDistributionPlot(y_sim_unknown,"blue",ann=true)
histogram!(y_sim_adult, bins=100, normalize=:pdf, alpha=0.1, linealpha=0.1, color="black")
histogram!(y_sim_kid, bins=100, normalize=:pdf, alpha=0.1, linealpha=0.1, color="red")
Plots.savefig(projDir*"/figs/PP_unknown_Slice.pdf")


##############
#### b-Version 2) not knowing that it is a child, set a fixed fraction
# using fraction 0.5
weight = 0.5
postBeingKid = weight
areKids = postBeingKid .>= rand(N)
idx = Int.(ceil.(rand(N).*N))   # choose a random index which will pick from the simulated posterior
zlogy_sim_unknown = μ_0[idx] .+ ϕ[idx] .* areKids .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
y_sim_unknown = exp.(zlogy_sim_unknown .* logStd .+ logMean)
makeDistributionPlot(y_sim_unknown,"blue",ann=true)
histogram!(y_sim_adult, bins=100, normalize=:pdf, alpha=0.1, linealpha=0.1, color="black")
histogram!(y_sim_kid, bins=100, normalize=:pdf, alpha=0.1, linealpha=0.1, color="red")
plot!(ann=(1200,0.0025,"Kids: $(weight*100) %"),grid=false)
Plots.savefig(projDir*"/figs/PP_unknown_fixed05_Slice.pdf")



#cur_colors = get_color_palette(:auto, plot_color(:white), 11)
# Limits for the figure
myXlims = (50,1200)
#myYlims = (0,0.006)

plt = StatsPlots.plot(size = (800, 1600),xlims=myXlims)

# Try various fixed ratios:
for (i, weight) in enumerate(0:.1:1)
    postBeingKid = weight
    areKids = postBeingKid .>= rand(N)
    idx = Int.(ceil.(rand(N).*N))   # choose a random index which will pick from the simulated posterior
    zlogy_sim_unknown = μ_0[idx] .+ ϕ[idx] .* areKids .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
    y_sim_unknown = exp.(zlogy_sim_unknown .* logStd .+ logMean)

    makeDistributionPlot!(y_sim_unknown,"blue",ann=false,offset=(i-1)*0.006,scale=1.0)
    makeDistributionPlot!(y_sim_adult,"black",ann=false,offset=(i-1)*0.006,scale=1.0)
    makeDistributionPlot!(y_sim_kid,"red",ann=false,offset=(i-1)*0.006,scale=1.0)
    plot!(ann=(myXlims[2]-200,(i-1)*0.006+0.003,"Kids: $(weight*100) %"),grid=false,yticks=false)
end
# Show plot
plt

Plots.savefig(projDir*"/figs/CompMixture_Slice.pdf")
