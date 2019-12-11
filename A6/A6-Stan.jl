#using Revise
using CmdStan
using Plots
using StatsPlots
using AxisArrays

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
405, 397, 383, 360, 387, 429, 358, 459, 371, 368, 452, 358, 371]

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

# Indicator for each individual j whether he/she is a child or not:
child_j = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1,
           0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0];

# Indicator for each observations i whether it comes from a child or not:
child_i = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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


#%% Functions ##################################################################
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

# A function to make the plots with HDIs
function makeDistributionPlot(X)
    #%% Mode, Mean, HDIs
    ω = mean(find_mode(X))
    μ_bar = mean(X);
    left,right = hdi(X);

    #%% Plotting
    # specify statsplot?
    histogram(X, bins=100, normalize=:pdf, label="MCMC", alpha=0.3, linealpha=0.1)   # Comes from StatsPlots now
    dPlot = density!(X,linewidth=3,label="density estimate")
    dCurve = filter(!isnan,dPlot.series_list[1].plotattributes[:y])
    dTopPoint = maximum(dCurve)
    plot!([(left,0),(left,dTopPoint/2)], linewidth=3, color="green", label="HDI",
           annotations = (left, dTopPoint/2, text("$(round(left,sigdigits=4))",:green,:bottom)))
    plot!([(right,0),(right,dTopPoint/2)],linewidth=3,color="green",label="HDI",
           annotations = (right, dTopPoint/2, text("$(round(right,sigdigits=4))",:green,:bottom)))
    plot!([(μ_bar,0),(μ_bar,dTopPoint)], linewidth=3, color="blue", label="mean",
           annotations = (μ_bar, 0, text("$(round(μ_bar,sigdigits=4))",:blue,:top)))
    plot!([(ω,0),(ω,dTopPoint)], linewidth=3, color="red", label="mode",
           annotations = (ω, dTopPoint, text("$(round(ω,sigdigits=4))",:red,:bottom)))

end

#%% Logarithmise and Transform Data ############################################
J = maximum(ind)    # number of individuals
I = length(y)       # number of observations

# Go to log-space:
logy = log.(y)

logMean = mean(logy)
logStd = std(logy)

# Mean-centre and scale each element:
zlogy = (logy .- logMean) ./ logStd;


################################################################################
############################# STANS MCMC #######################################

#%% Stan Setup #################################################################

# If run from command line externally:
#projDir = dirname(@__FILE__)

# If run from Jupyter/Hydrogen, maybe change to suit you:
projDir= "/home/johhub/Desktop/ABDA/A6"
#projDir= "/lhome/johhub/Desktop/ABDA/A5"
tmpDir = projDir*"/tmp"

noOfChains = 4
N = 10^5 / noOfChains   # more than 10^6 samples make the histograms thin
N = convert(Int64, N)
keepchains = false
burnIn = 10^4

modelString = "
data {
    int<lower=1> J;         // number of individuals for which data was observed
    int<lower=1> I;         // observation number in the total vector
    int<lower=1> ID[I];     // individual ID in the total vector
    int<lower=0> ISKID[J];  // identifief for whether individual is a kid
    real zlogy[I];           // reaction time measurements in the total vector
}
parameters {
    real theta[J];          // mean for each individual
    real<lower=0.000001> sigma;    // same std for all individuals
    real mu;                // mean for the group
    real<lower=0.000001> tau;      // std of the group (only on group = all individuals)
    real phi;
}
transformed parameters { // no transformed variables to use
// how to transform here?
}
model {
    for (i in 1:I)
        zlogy[i] ~ normal(theta[ID[i]],sigma);
    for (j in 1:J)
        theta[j] ~ normal(mu + phi*ISKID[j],tau);
    // no prior is equivalent to a uniform prior
    // however this is an <<improper>> prior and can lead to problems
}
generated quantities {
    // real zlogy_pred;
    // real theta_pred;
    // theta_pred = normal_rng(mu,tau);
    // zlogy_pred = normal_rng(theta_pred,sigma);
}
";


#%%### 2) Data for the Stan model, note that variable names must
# correspond to defined model in stan-code
observedData = Dict("I" => I,
                    "J" => J,
                    "zlogy" => zlogy,
                    "ID" => ind,
                    "ISKID" => child_j);

### 3) Chain specs
myModel = Stanmodel(
                Sample(save_warmup=false,
                       num_warmup=burnIn,
                       num_samples=N,
                       thin=1),   # thin: Period between saved samples
                name = "reactionTime-A6",
                model = modelString,
                printsummary = false,
                tmpdir = tmpDir,
                nchains = noOfChains);    # number of chains (4 default)

#%%### 4) Run CmdStan:
# rc:     return code (0 if all is fine)
# chn:    chain results
# cnames: vector of variable names
rc, chn, cnames = stan(myModel,
                       observedData,
                       projDir,
                       diagnostics = false,
                       CmdStanDir = CMDSTAN_HOME);


#%% Transform back (undo mean-centering and scaling and go to non-log space)
#
# we assumed mean-centered and scaled log-normal distributions for theta
# by our "tricks" above, i.e. μ_y, σ_y fixed
#                    z_logy   ~  Normal(θ, σ^2)
#       (log(y) - m_y) / σ_y  ~  Normal(θ, σ^2)
#       (log(y) - m_y) / σ_y  =  θ + σ * x
# where y is the observation and x the predictor,
# and x is Normal(0,1), then
#       log(y) = (θ + σ * x) * σ_y + m_y
#       log(y) = θ*σ_y + σ*σ_y*x  + m_y
#       log(y) = θ*σ_y + m_y + σ*σ_y*x
#
#  say θ_trans = θ*σ_y + m_y
#      σ_trans = σ*σ_y
#
#       log(y) = θ_trans + σ_trans * x
#       log(y) ~ Normal(θ_trans, σ_trans^2)
#           y  = exp(θ_trans + σ_trans * x)
#           y  ~ LogNormal(θ_trans, σ_trans^2)
#         E[y] = exp(θ_trans + σ_trans^2 / 2)   # from Wiki
#        which is what θ, the individual mean, expresses
#
# The same is true for the prior distribution yielding θ:
# (μ and τ are hyperparameters):
#                         θ ~ Normal(μ_0 + ϕ,τ^2)   # for kids
#                         θ ~ Normal(μ_0 + 0,τ^2)   # for adults
#     (θ_trans - m_y) / σ_y = (μ_0 + ϕ[0,1]) + τ * ξ
#   where ξ ~ Normal(0,1)
#
#                   θ_trans = (μ_0 + ϕ[0,1])*σ_y + m_y + τ*σ_y*ξ
#   say μ_0_trans =  μ_0*σ_y + m_y
#       μ_ϕ_trans = (μ_0 + ϕ)*σ_y + m_y
#         τ_trans = τ*σ_y
#                   θ_trans = μ_trans + τ_trans*ξ
#                   θ_trans ~ Normal(μ_trans, τ_trans^2)
#
#                    log(y) = μ_trans + τ_trans*ξ + σ_trans*x
# ξ and x are both ~ N(0,1)
# a sum of two distributions, i.e. here N(0,τ_trans^2) + N(0,σ_trans^2), are a
# convolution of the two, thus their sum is: N(0,τ_trans^2 + σ_trans^2)
# assuming independence, and therefore
#                    log(y) = μ_trans + sqrt(τ_trans^2 + σ_trans^2)*x
#                    log(y) ~ Normal(μ_trans, (σ_trans^2 + τ_trans^2))
#                      E[y] = exp(μ_trans + (σ_trans^2 + τ_trans^2) / 2)
#                      E[y] = exp(μ_trans + σ_trans^2 / 2 + τ_trans^2 / 2)
#                 which is what μ, the group mean, expresses

#%% check the names and positions of the vars (might change with naming)
chn.value[1,:,1]

#%%
# Originals
# indices in chains:
# 5 = mu
# 7 = sigma
# 9 = tau
# 10 = theta1
# 43 = theta34
# 44 = theta_pred
# 46 = logy_pred

# Access the axis array like this:
ϕ = 1.0 * chn.value[Axis{:var}("phi")][:]   # all chains in one sausage
θ = Array{Float64,2}(undef,N*noOfChains,J)
for j in 1:J
    θ[:,j] = 1.0 * chn.value[Axis{:var}("theta.$j")][:]   # all chains in one sausage
end
μ = 1.0 * chn.value[Axis{:var}("mu")][:]
σ = 1.0 * chn.value[Axis{:var}("sigma")][:]
τ = 1.0 * chn.value[Axis{:var}("tau")][:]

makeDistributionPlot(ϕ)

Plots.savefig("/home/johhub/Desktop/ABDA/A6/test-Stan.pdf")

makeDistributionPlot(μ)

#logy_pred = 1.0 * chn.value[:,46,1]


#%%
# Un-scale and un-mean-centre:
θ_trans = θ .* logStd .+ logMean
μ_0_trans = (μ .+ 0) .* logStd .+ logMean
# kids only:
μ_ϕ_trans = (μ .+ ϕ) .* logStd .+ logMean

σ_trans = σ .* logStd
τ_trans = τ .* logStd

#logy_trans = logy_pred .* logStd .+ logMean
#logy_trans = logy_pred
ϕ_trans = ϕ .*  logStd .+ logMean

# Get into non-log space:
θ_trans_unLog = exp.(θ_trans .+ 0.5 .* repeat(σ_trans,1,J).^2);

μ_0_trans_unLog = exp.(μ_0_trans .+ 0.5 .* σ_trans.^2 .+ 0.5 .* τ_trans.^2);
μ_ϕ_trans_unLog = exp.(μ_ϕ_trans .+ 0.5 .* σ_trans.^2 .+ 0.5 .* τ_trans.^2);

# need to transform μ+ϕ together? instead, then subtract the mu
#ϕ_trans_unLog = exp.(ϕ_trans .+ 0.5 .* σ_trans.^2 .+ 0.5 .* τ_trans.^2);
#ϕ_trans_unLog = μ_ϕ_trans_unLog - μ_0_trans_unLog;
ϕ_trans_unLog = logStd.*ϕ


#τ_trans_unLog = sqrt.((exp.(τ_trans.^2) .- 1.0) .* (2.0 .* μ_0_trans .+ τ_trans.^2))
τ_trans_unLog = sqrt.((exp.(τ_trans.^2) .- 1.0) .* (2.0 .* μ_ϕ_trans .+ τ_trans.^2))
#logy_trans_unLog = exp.(logy_trans);
makeDistributionPlot(μ_trans_unLog)
makeDistributionPlot(θ_trans_unLog[:,1])

################################################################################
############################# TASKS ############################################

######################## Task 1 ####################################
makeDistributionPlot(ϕ)
makeDistributionPlot(ϕ_trans_unLog)
# Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A1-Dude-Stan.pdf")


# ######################## Task 2 #############
makeDistributionPlot(τ)
makeDistributionPlot(τ_trans_unLog)
# Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Group-Stan.pdf")
#
#
# ########################  Task A-2-a-ii (new individual prediction)  ###########
# samplePts = logy_trans_unLog
# ω = mean(find_mode(samplePts))
# μ_bar = mean(samplePts);
# med = median(samplePts);
# left,right = hdi(samplePts);
#
# histogram(samplePts, bins=100, normalize=:pdf, label="Simulation", alpha=0.3, linealpha=0.1)   # Comes from StatsPlots now
# dPlot = density!(samplePts,linewidth=3,label="density estimate")
# dCurve = filter(!isnan,dPlot.series_list[1].plotattributes[:y])
# dTopPoint = maximum(dCurve)
# plot!([(left,0),(left,dTopPoint/2)], linewidth=3, color="green", label="HDI",
#        annotations = (left, dTopPoint/2, text("$(Int(round(left)))",:green,:bottom)))
# plot!([(right,0),(right,dTopPoint/2)],linewidth=3,color="green",label="HDI",
#        annotations = (right, dTopPoint/2, text("$(Int(round(right)))",:green,:bottom)))
# plot!([(μ_bar,0),(μ_bar,dTopPoint)], linewidth=3, color="blue", label="mean",
#        annotations = (μ_bar, 0, text("$(Int(round(μ_bar)))",:blue,:top)))
# plot!([(ω,0),(ω,dTopPoint)], linewidth=3, color="red", label="mode",
#        annotations = (ω, dTopPoint, text("$(Int(round(ω)))",:red,:bottom)))
# plot!([(med,0),(med,dTopPoint)], linewidth=3, color="purple", label="median",
#        annotations = (med, 0, text("$(Int(round(med)))",:purple,:bottom)))
#
# Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Pred-Stan.pdf")
#
#
# ########################  Task A-2-b (compare to website statistics)  ##########
# #%% Taken from the homepage:
# median_web = 273.0
# mean_web = 284.0
# println("Δ Median = ",median_web - med)  # using median from previous task
# println("Δ Mean = ",mean_web - μ_bar)    # using mean from previous task
#
#
# ########################  Task A-3 #############################################
# # Compare hierarchical theta to individual theta using sample means
#
# #%% Logarithmic sample means
# θ_means_log = zeros(J)
# for j in 1:J
#     θ_means_log[j] = mean(logy[ind .== j])
# end
#
# #%% Plot into one figure
# # Get ordered indices:
# sortIdx = sortperm(θ_means_log)
# # Limits for the figure
# myXlims = (minimum(θ_trans),maximum(θ_trans))
# # Initialise the subplots
# StatsPlots.plot(layout=(J, 1),size = (1000, 1500))
# # Plot each theta:
# for i in 1:(J-1)
#     j = sortIdx[i]
#     # Sampled thetas and their mean:
#     histogram!(θ_trans[:,j], bins=100, normalize=:pdf, legend=false, alpha=0.3, linealpha=0.0,
#             ann=(myXlims[1]+.05,4,"ind $j:"),ticks=nothing, yaxis=false, subplot=i, xlims=myXlims)
#     vline!([mean(θ_trans[:,j])],linewidth=3, color="black", subplot=i, legend=false)
#
#     # The (log) sample means of the initial data:
#     vline!([θ_means_log[j]],linewidth=3, color="red", subplot=i, legend=false)
# end
#
# # The last one separately so I can see it in Hydrogen:
# j = sortIdx[J]
# histogram!(θ_trans[:,j], bins=100, normalize=:pdf, legend=false, alpha=0.3, linealpha=0.0,
#             ann=(myXlims[1]+.05,10,"ind $J:"),ticks=nothing, yaxis=false, subplot=J, xlims=myXlims)
# vline!([mean(θ_trans[:,j])],linewidth=3, color="black", subplot=J, legend=false)
# vline!([θ_means_log[j]],linewidth=3, color="red", subplot=J, legend=false)
#
# Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Comp-MLE-Slice-All.pdf")
#
#
# # --------------- Old Code (please ignore) -------------------------
# # Get a color map:
# #curColor = get_color_palette(:auto, plot_color(:white), J)
# #for n = 1:Int(ceil(J/5))
# #    plot()
# #    for j in ((n-1)*5+1):min((n*5),J)
# #        histogram!(θ_trans[:,j], bins=100, normalize=:pdf,legend=false,alpha=0.2, linealpha=0.0, color=curColor[j])
# #        vline!([θ_means_log[j]], linewidth=3, color=curColor[j])
# #    end
# #    Plots.savefig("/home/johhub/Desktop/ABDA/A5/figs/A2-Comp-MLE-Stan-$n.pdf")
# #end
