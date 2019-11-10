using Revise
using CmdStan
using Plots
using StatsPlots
using Random

################################################################################
############################# PREPARATION ######################################

#%% Input Data #################################################################
y = [607, 583, 521, 494, 369, 782, 570, 678, 467, 620, 425, 395, 346, 361, 310,
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

ind = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 5, 5, 5, 5,
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


#%% Functions ##################################################################
function mySliceSampler(pdf_log, x_0, w, m, N = 2000, burnIn = 1000)
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

# A recursive application of the HDI algorithm to find the mode:
function find_mode(sample_pts;runs=2)
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

function log_prior_flat(X::Float64)
    return 0.0  # or any other constant value, for the posterior it doesn't matter
end


function log_prior_flat_positiveOnly(X::Float64)
    if X > 0
        return 0.0  # or any other constant value, for the posterior it doesn't matter
    else
        return -Inf
    end
    # Alternative of Jesper:
    #return log(tau > 0)
end


function log_normal(θ::Array{Float64},μ::Float64,τ::Float64)
    # Return the joint pdf of an array of thetas, i.e. the sum in log-space

    # We can leave out the 1/sqrt(2π) since it does not change sampling from the posterior
    return sum(-log(τ) .- 0.5 .* ((θ .- μ) ./ τ).^2)
end


function log_lklhd_fct(jointParams::Array{Float64},zlogy::Array{Float64,1},ind::Array{Int64,1})
    # return the joint lklhd of the given input data
    J = length(jointParams) - 3
    θ = jointParams[1:J]
    μ, σ, τ  = jointParams[(J+1):end]

    jointLklhd = 0.0

    if (σ > 0)
        for i in 1:length(zlogy)
            jointLklhd += -log(σ) - 0.5 * ((zlogy[i] - θ[ind[i]]) / σ)^2
        end
    else
        jointLklhd = -Inf
    end

    return jointLklhd
end


function log_priors(jointParams::Array{Float64})
    # The joint log-pdf of all the priors in the model
    J = length(jointParams) - 3
    θ = jointParams[1:J]
    μ, σ, τ  = jointParams[(J+1):end]

    if (σ > 0) && (τ > 0)
        return log_normal(θ,μ,τ) + log_prior_flat(μ) + log_prior_flat_positiveOnly(σ) + log_prior_flat_positiveOnly(τ)
    else
        return -Inf
    end
end


#%% Logarithmise and Transform Data ############################################
J = maximum(ind)    # number of individuals
I = length(y)       # number of observations

# Go to log-space:
logy = log.(y)

logMean = mean(logy)
logStd = std(logy)

# Mean-centre and scale each element:
zlogy = (logy .- logMean) ./ logStd


#%% Concretise the functions by including the data #############################
log_lklhd(jointParams::Array{Float64}) = log_lklhd_fct(jointParams::Array{Float64}, zlogy, ind)
log_posterior(jointParams::Array{Float64}) = log_lklhd(jointParams::Array{Float64}) + log_priors(jointParams::Array{Float64})


################################################################################
############################# MCMC #############################################

#%% Chain Setup
N = 10^6
burnIn = 10^4
x_start = ones(Float64,J+3) * 0.1
w = ones(Float64,J+3) * 0.1
m = 100

#%% MCMC run
@time chain, pdf = mySliceSampler(log_posterior,x_start,w,m,N,burnIn);

#%% Transform back (undo mean-centering and scaling and go to non-log space)
#
# we assumed mean-centered and scaled log-normal distributions for theta
# by our "tricks" above, i.e. μ_y, σ_y fixed
#                    z_logy   ~  Normal(θ, σ^2)
#       (log(y) - m_y) / σ_y  ~  Normal(θ, σ^2)
#       (log(y) - m_y) / σ_y  =  θ + σ * x
# where y is the observation and x the predictor, then
#       log(y) = (θ + σ * x) * σ_y + m_y
#       log(y) = θ*σ_y + σ*σ_y*x  + m_y
#       log(y) = θ*σ_y + m_y + σ*σ_y*x
#       log(y) = θ_trans + σ_trans * x
#       log(y) ~ Normal(θ_trans, σ_trans^2)
#           y  = exp(θ_trans + σ_trans * x)
#           y  ~ LogNormal(θ_trans, σ_trans^2)
#         E[y] = exp(θ_trans + σ_trans^2 / 2)
#
# The same is true for the hyperdistribution yielding theta:
#                         θ ~ Normal(μ,τ^2)
#     (θ_trans - m_y) / σ_y = μ + τ * xi
#                   θ_trans = μ*σ_y + m_y + τ*σ_y*xi
#                   θ_trans = μ_trans + τ_trans*xi
#                   θ_trans ~ Normal(μ_trans, τ_trans^2)
#
#                    log(y) ~ μ_trans + τ_trans*xi + σ_trans*x
# xi and x are both ~ N(0,1)
# a sum of two distributions, i.e. here N(0,τ_trans^2) + N(0,σ_trans^2), are a
# convolution of the two, thus their sum is: N(0,τ_trans^2 + σ_trans^2)
# assuming independence, and therefore
#                    log(y) ~ μ_trans + sqrt(τ_trans^2 + σ_trans^2)*x
#                    log(y) ~ Normal(μ_trans, (σ_trans^2 + τ_trans^2))
#                      E[y] = exp(μ_trans + (σ_trans^2 + τ_trans^2) / 2)
#                      E[y] = exp(μ_trans + σ_trans^2 / 2 + τ_trans^2 / 2)

# Originals
θ = chain[:,1:J]
μ = chain[:,(J+1)]
σ = chain[:,(J+2)]
τ = chain[:,(J+3)]

# Un-scale and un-mean-centre:
θ_trans = θ .* logStd .+ logMean
μ_trans = μ .* logStd .+ logMean
σ_trans = σ .* logStd
τ_trans = τ .* logStd

# Get into non-log space:
θ_trans_unLog = exp.(θ_trans .+ 0.5 .* repeat(σ_trans,1,J).^2)
μ_trans_unLog = exp.(μ_trans .+ 0.5 .* σ_trans.^2 .+ 0.5 .* τ_trans.^2)


################################################################################
############################# TASKS ############################################

######################## Task A-1 (the Dude)####################################
#%% Mode, Mean, HDIs
dude = θ_trans_unLog[:,4]
ω = mean(find_mode(dude))
μ_bar = mean(dude);
left,right = hdi(dude);

#%% Plotting
histogram(dude, bins=100, normalize=:pdf, label="MCMC", alpha=0.3, linealpha=0.1)   # Comes from StatsPlots now
dPlot = density!(dude,linewidth=3,label="density estimate")
dCurve = filter(!isnan,dPlot.series_list[1].plotattributes[:y])
dTopPoint = maximum(dCurve)
plot!([(left,0),(left,dTopPoint/2)], linewidth=3, color="green", label="HDI",
       annotations = (left, dTopPoint/2, text("$(Int(round(left)))",:green,:bottom)))
plot!([(right,0),(right,dTopPoint/2)],linewidth=3,color="green",label="HDI",
       annotations = (right, dTopPoint/2, text("$(Int(round(right)))",:green,:bottom)))
plot!([(μ_bar,0),(μ_bar,dTopPoint)], linewidth=3, color="blue", label="mean",
       annotations = (μ_bar, 0, text("$(Int(round(μ_bar)))",:blue,:top)))
plot!([(ω,0),(ω,dTopPoint)], linewidth=3, color="red", label="mode",
       annotations = (ω, dTopPoint, text("$(Int(round(ω)))",:red,:bottom)))

Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A1-Dude-Slice.pdf")


######################## Task A-2-a-i (expectation new individual) #############
#%% Mode, Mean, HDIs
samplePts = μ_trans_unLog
ω = mean(find_mode(samplePts))
μ_bar = mean(samplePts);
left,right = hdi(samplePts);

#%% Plotting
histogram(samplePts, bins=100, normalize=:pdf, label="MCMC", alpha=0.3, linealpha=0.1)   # Comes from StatsPlots now
dPlot = density!(samplePts,linewidth=3,label="density estimate")
dCurve = filter(!isnan,dPlot.series_list[1].plotattributes[:y])
dTopPoint = maximum(dCurve)
plot!([(left,0),(left,dTopPoint/2)], linewidth=3, color="green", label="HDI",
       annotations = (left, dTopPoint/2, text("$(Int(round(left)))",:green,:bottom)))
plot!([(right,0),(right,dTopPoint/2)],linewidth=3,color="green",label="HDI",
       annotations = (right, dTopPoint/2, text("$(Int(round(right)))",:green,:bottom)))
plot!([(μ_bar,0),(μ_bar,dTopPoint)], linewidth=3, color="blue", label="mean",
       annotations = (μ_bar, 0, text("$(Int(round(μ_bar)))",:blue,:top)))
plot!([(ω,0),(ω,dTopPoint)], linewidth=3, color="red", label="mode",
       annotations = (ω, dTopPoint, text("$(Int(round(ω)))",:red,:bottom)))

Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Group-Slice.pdf")


########################  Task A-2-a-ii (new individual prediction)  ############
#%% Predict a reaction time for a single measurement (y)
zlogy_sim = zeros(N)
for i in 1:N
    # 1) Pick a posterior sample from mu, tau, and sigma
    idx = Int(ceil(rand()*N))   # choose a random index which will pick from the simulated posterior
    μ_pick = μ[idx]
    τ_pick = τ[idx]
    σ_pick = σ[idx]
    # 2) Simulate a new theta given these samples, i.e. theta ~ N(μ,τ)
    θ_sim  = μ_pick + randn()*τ_pick
    # 3) Simulate a reaction time measurement (y) given the picked theta and sigma, i.e. zlogy ~ N(θ,σ)
    zlogy_sim[i] = θ_sim + randn()*σ_pick
end

# 4) Calculate y from zlogy (undo mean-centering and scaling, and go to non-log space)
y_sim = exp.(zlogy_sim .* logStd .+ logMean)
samplePts = y_sim

ω = mean(find_mode(samplePts))
μ_bar = mean(samplePts);
med = median(samplePts);
left,right = hdi(samplePts);

histogram(samplePts, bins=100, normalize=:pdf, label="Simulation", alpha=0.3, linealpha=0.1)   # Comes from StatsPlots now
dPlot = density!(samplePts,linewidth=3,label="density estimate")
dCurve = filter(!isnan,dPlot.series_list[1].plotattributes[:y])
dTopPoint = maximum(dCurve)
plot!([(left,0),(left,dTopPoint/2)], linewidth=3, color="green", label="HDI",
       annotations = (left, dTopPoint/2, text("$(Int(round(left)))",:green,:bottom)))
plot!([(right,0),(right,dTopPoint/2)],linewidth=3,color="green",label="HDI",
       annotations = (right, dTopPoint/2, text("$(Int(round(right)))",:green,:bottom)))
plot!([(μ_bar,0),(μ_bar,dTopPoint)], linewidth=3, color="blue", label="mean",
       annotations = (μ_bar, 0, text("$(Int(round(μ_bar)))",:blue,:top)))
plot!([(ω,0),(ω,dTopPoint)], linewidth=3, color="red", label="mode",
       annotations = (ω, dTopPoint, text("$(Int(round(ω)))",:red,:bottom)))
plot!([(med,0),(med,dTopPoint)], linewidth=3, color="purple", label="median",
       annotations = (med, 0, text("$(Int(round(med)))",:purple,:bottom)))

Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Pred-Slice.pdf")


########################  Task A-2-b (compare to website statistics)  ##########
#%% Taken from the homepage:
median_web = 273.0
mean_web = 284.0
println("Δ Median = ",median_web - med)
println("Δ Mean = ",mean_web - μ_bar)


########################  Task A-3 #############################################
#%% Compare hierarchical theta to individual theta using sample means
histogram(θ_trans_unLog[:,1], bins=100, normalize=:pdf,alpha=0.2, linealpha=0.0)
n = 10
for i in 2:34
    histogram!(θ_trans_unLog[:,i], bins=100, normalize=:pdf, alpha=0.2, linealpha=0.0)
end
histogram!(θ_trans_unLog[:,n+1], bins=100, normalize=:pdf, alpha=0.2, linealpha=0.0)
histogram!(μ_trans_unLog, bins=100, normalize=:pdf, alpha=1, linealpha=0.1)

#%% New model without hierarchy
# Logarithmic sample means
θ_means_log = zeros(J)
#for i in 1:I
#    n = length(findall(x -> x==ind[i],ind)) # finding the number of data points for each individual
#    θ_means_log[ind[i]] += logy[i] / n  # sum up the weighted results
#end
# or use this alternative:
for j in 1:J
    θ_means_log[j] = mean(logy[ind .== j])
end

θ_means_log
#%%
# theta ~ N(θ_means_log,sigma)

#histogram(θ_trans[:,1], bins=100, normalize=:pdf,label="theta from MCMC",alpha=0.2, linealpha=0.0)
#plot([(θ_means_log[1],0),(θ_means_log[1],.5)], linewidth=3, color="red", label="sample mean",
#       annotations = (med, 0, text("$(Int(round(θ_means_log[1])))",:red,:bottom)))
#plot()
#for i in 1:12
#    histogram!(θ_trans[:,i], bins=100, normalize=:pdf, legend=false, alpha=0.3, linealpha=0.0,
#            inset = (1, bbox(0.0,0.1*i,1,0.1,:bottom,:left)),  # x,y, width, height, origin
#            ticks=nothing, subplot=1+i, bg_inside=nothing)
#vline!([θ_means_log[i]],linewidth=3, color="red", subplot=i+1, legend=false)
#end

# Get a color map:
curColor = get_color_palette(:auto, plot_color(:white), J)
for n = 1:Int(ceil(J/5))
    plot()
    for j in ((n-1)*5+1):min((n*5),J)
        histogram!(θ_trans[:,j], bins=100, normalize=:pdf,legend=false,alpha=0.2, linealpha=0.0, color=curColor[j])
        vline!([θ_means_log[j]], linewidth=3, color=curColor[j])
    end

Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Comp-MLE-Stan-$n.pdf")

end

#%%

histogram(θ_trans, layout=(5,1), bins=100, normalize=:pdf,legend=false,alpha=0.2, linealpha=0.0)
vline!(θ_means_log, layout=(5,1), linewidth=3)
