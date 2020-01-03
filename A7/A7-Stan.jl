using CmdStan
using Plots
using StatsPlots
using AxisArrays
using StatsBase
using Random
using Distributions

##### Some local setups
# If run from command line externally:
#projDir = dirname(@__FILE__)

# If run from Jupyter/Hydrogen, maybe change to suit you:
projDir= "/lhome/johhub/Desktop/ABDA/A7"
tmpDir = projDir*"/tmp"

#####
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

# Indicator for each observation i whether it comes from a child or not:
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

# Attempt numbers of each observation i:
x = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 1, 2, 3,
4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3,
4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
20, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2,
3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 2, 3, 4,
5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5,
6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4, 5, 6, 7, 8,
9, 10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1,
2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 1, 2, 3, 4,
5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7,
8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 1, 2, 3, 4,
5, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18];

# Alternatively, have it all in one array:
x_alt = Array{Int64,1}[
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5, 6, 7, 8, 9],
[1], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16],
[1, 2, 3, 4, 5, 6, 7, 8, 9],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
[1, 2, 3, 4, 5, 6],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5, 6, 7, 8, 9],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
[1, 2, 3, 4, 5, 6],
[1, 2, 3, 4, 5],
[1, 2, 3, 4, 5],
[1],
[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]];


#####
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
function makeDistributionPlot(X, color="blue"; ann=true, offset=0.0, scale=1.0)
    #%% Mode, Mean, HDIs
    ω = mean(find_mode(X))
    μ_bar = mean(X);
    left,right = hdi(X);

    h = fit(Histogram, X, nbins=200)
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

    h = fit(Histogram, X, nbins=200)
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

#%% Logarithmise and Transform Data ############################################
J = maximum(ind)    # number of individuals
I = length(y)       # number of observations

# Go to log-space:
logy = log.(y)

logMean = mean(logy)
logStd = std(logy)

global trainMean = mean(x)
global trainStd = std(x)

# Mean-centre and scale each element:
zlogy = (logy .- logMean) ./ logStd;

# Standardise training data:
zx = (x .- trainMean) ./ trainStd;


################################################################################
############################# STANS MCMC #######################################

#%% Stan Setup #################################################################
noOfChains = 4
N = 10^5                # total number of samples
N_chain = convert(Int64, N / noOfChains)   # samples per chain
keepchains = false
burnIn = 10^3           # burn-in per chain

modelString = "
data {
    int<lower=1> J;         // number of individuals for which data was observed
    int<lower=1> I;         // observation number in the total vector
    int<lower=1> ID[I];     // individual ID in the total vector
    int<lower=0> K[J];      // input: identifier for whether individual is a kid
    real zx[I];             // (transformed) input: attempt number of the observation
    real zlogy[I];          // (transformed) output: reaction time
}
parameters {
//    real theta0[J];          // intercept for each individual
//    real theta1[J];          // slope for each individual
    real eta0[J];          // helper for reparametrisation
    real eta1[J];          // helper for reparametrisation
    real<lower=0.000001> sigma;    // variation around regression line (same for all)
    real mu0;               // mean for the intercept at group level
    real phi0;              // additive term for the intercept if being a kid
    real mu1;               // mean for the slope at the group level
    real phi1;              // additive term for the slope if being a kid
    real tau0;     // std of the group intercept (same for all)
    real tau1;   // std of the group slope (same for all)
}
transformed parameters {
    real theta0[J];          // intercept for each individual
    real theta1[J];          // slope for each individual
    for (j in 1:J) {
        theta0[j] = mu0 + phi0*K[j] + tau0 * eta0[j];
        theta1[j] = mu1 + phi1*K[j] + tau1 * eta1[j];
    }
}
model {
    tau0 ~ uniform(0,10000);
    tau1 ~ uniform(0,10000);
    for (i in 1:I)
        zlogy[i] ~ normal(theta0[ID[i]] + theta1[ID[i]]*zx[i], sigma);
    for (j in 1:J) {
        eta0[j] ~ normal(0, 1);
        eta1[j] ~ normal(0, 1);
    }

    //for (j in 1:J) {
    //    theta0[j] ~ normal(mu0 + phi0*K[j], tau0);
    //    theta1[j] ~ normal(mu1 + phi1*K[j], tau1);
    //}
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
                    "K" => child_j,
                    "zx" => zx);

### 3) Chain specs
myModel = Stanmodel(
                Sample(save_warmup = false,
                       num_warmup = burnIn,
                       num_samples = N_chain,
                       thin = 1),   # thin: Period between saved samples
                name = "reactionTime-A7",
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

# takes for N=10^4 around 30 seconds.
# N=10^5 around 1.5 minutes


################################################################################
############################# RESULTS ##########################################

##############
# Acess results and Transform back (undo mean-centering and scaling and go to non-log space)
# See derivation in PDF file instead

# Access the axis array by names:
θ_0 = Array{Float64,2}(undef,N,J)
θ_1 = Array{Float64,2}(undef,N,J)
for j in 1:J
    θ_0[:,j] = 1.0 * chn.value[Axis{:var}("theta0.$j")][:]
    θ_1[:,j] = 1.0 * chn.value[Axis{:var}("theta1.$j")][:]
end
μ_0 = 1.0 * chn.value[Axis{:var}("mu0")][:]
μ_1 = 1.0 * chn.value[Axis{:var}("mu1")][:]
ϕ_0 = 1.0 * chn.value[Axis{:var}("phi0")][:]   # all chains in one sausage
ϕ_1 = 1.0 * chn.value[Axis{:var}("phi1")][:]
σ   = 1.0 * chn.value[Axis{:var}("sigma")][:]
τ_0 = 1.0 * chn.value[Axis{:var}("tau0")][:]
τ_1 = 1.0 * chn.value[Axis{:var}("tau1")][:]


##### RE-READ RESULTS IF CHAIN HAS ALREADY RUN
# only for re-use purposes, please ignore!
using CSV

θ_0 = Array{Float64,2}(undef,N,J)
θ_1 = Array{Float64,2}(undef,N,J)
μ_0 = Array{Float64,1}(undef,N)
μ_1 = Array{Float64,1}(undef,N)
ϕ_0 = Array{Float64,1}(undef,N)
ϕ_1 = Array{Float64,1}(undef,N)
σ   = Array{Float64,1}(undef,N)
τ_0 = Array{Float64,1}(undef,N)
τ_1 = Array{Float64,1}(undef,N)

for i in 1:noOfChains
    # subindices
    from = N_chain*(i-1) + 1
    to = N_chain*i

    # read into a DataFrame:
    chain_csv = CSV.read(tmpDir*"/reactionTime-A7_samples_$i.csv"; comment="#", normalizenames=true)

    θ_0_df = chain_csv[:,r"theta0"]    # regex
    θ_1_df = chain_csv[:,r"theta1"]
    for j in 1:J
        θ_0[from:to,j] = θ_0_df[:,j]  # need to bring on array form, from data frame
        θ_1[from:to,j] = θ_1_df[:,j]  # need to bring on array form, from data frame
    end

    μ_0[from:to,:] = chain_csv[:, r"mu0"][:,1]
    μ_1[from:to,:] = chain_csv[:, r"mu1"][:,1]
    ϕ_0[from:to,:] = chain_csv[:, r"phi0"][:,1]
    ϕ_1[from:to,:] = chain_csv[:, r"phi1"][:,1]
    σ[from:to,:]   = chain_csv[:, r"sigma"][:,1]
    τ_0[from:to,:] = chain_csv[:, r"tau0"][:,1]
    τ_1[from:to,:] = chain_csv[:, r"tau1"][:,1]
end
######

##### RE-READ RESULTS FROM OLD ASSIGNMENT
# using CSV
#
# A5_import = CSV.read(projDir*"/export_tau_A5.csv"; header=false)
# τ_A5 = A5_import[:,1]
# τ_A5_unscaled = A5_import[:,2]
#
# A5_import = CSV.read(projDir*"/export_mu_A5.csv"; header=false)
# μ_A5 = A5_import[:,1]
# μ_A5_unscaled = A5_import[:,2]
#####

#####
# Un-standardise:
# Individual level:
global θ_0_unscaled = logStd.*θ_0 .- logStd.*θ_1.*trainMean./trainStd .+ logMean
global θ_1_unscaled = logStd.*θ_1./trainStd

global σ_unscaled = σ .* logStd

# Group level (maybe not needed)
global μ_0_unscaled = logStd.*μ_0 .- logStd.*μ_1.*trainMean./trainStd .+ logMean
global ϕ_0_unscaled = logStd.*(ϕ_0 .- ϕ_1.*trainMean./trainStd)

global μ_1_unscaled = logStd.*μ_1./trainStd
global ϕ_1_unscaled = logStd.*ϕ_1./trainStd

global τ_0_unscaled = logStd.* sqrt.(τ_0.^2 .+ τ_1.^2 .* trainMean.^2 ./ trainStd.^2)
global τ_tilde_1_unscaled = sqrt.(2 .* logStd.^2 .* τ_1.^2 .* trainMean ./ trainStd.^2)
global τ_1_unscaled = logStd .* τ_1 .* trainMean ./ trainStd

#####
# Get into non-log space (generates function of input):
# Individual level:
function E_y_ind(x,j)
    # transform to zx:
    #zx = (x - trainMean) / trainStd
    return exp.(θ_0_unscaled[:,j] .+
                θ_1_unscaled[:,j] .* x .+
                0.5 .* σ_unscaled.^2);
end

function E_y_ind2(x,j)
    # transform to zx:
    #zx = (x - trainMean) / trainStd
    return exp.(θ_0[:,j] .+
                θ_1[:,j] .* x .+
                0.5 .* σ.^2);
end


#ϕ_0_unscaled_unLog = exp.(ϕ_unscaled)

#μ_0_trans_unLog = exp.(μ_0_unscaled .+
#                  0.5 .* σ_unscaled.^2 .+
#                  0.5 .* τ_0_unscaled.^2 .+
#                  0.5 .* τ_1_unscaled.^2 * zx);
#μ_ϕ_trans_unLog = exp.(μ_0_unscaled .+ ϕ_unscaled .+ 0.5 .* σ_unscaled.^2 .+ 0.5 .* τ_unscaled.^2);  # kids


################################################################################
############################# TASKS ############################################

makeDistributionPlot(ϕ_1,"blue")

makeDistributionPlot(τ_1,"blue")

makeDistributionPlot(exp.(θ_1_unscaled[:,4]),"blue")

#makeDistributionPlot(E_y_ind2(1,1),"blue")
#makeDistributionPlot!(E_y_ind2(5,1),"red")

makeDistributionPlot(E_y_ind2(1,3),"blue")
makeDistributionPlot!(E_y_ind2(5,3),"red")
plot!(xlims=[200,800])

attempts = range(0,22,length=100)
plot()
out = mean.(E_y_ind.(attempts,4))
plot!(attempts,out)

######################## Task 1 ####################################
# Effect of being a kid:
makeDistributionPlot(ϕ, "orange")
Plots.savefig(projDir*"/figs/phi_Stan.pdf")
makeDistributionPlot(ϕ_unscaled, "orange")
Plots.savefig(projDir*"/figs/phi_unscaled_Stan.pdf")

makeDistributionPlot(ϕ_unscaled_unLog, "orange")
Plots.savefig(projDir*"/figs/phi_unscaled_nonlog_Stan.pdf")


######################## Task 2 #############
makeDistributionPlot(τ,"blue")
makeDistributionPlot!(τ_A5,"purple")
Plots.savefig(projDir*"/figs/tau_Stan.pdf")

makeDistributionPlot(τ_unscaled, "blue")
makeDistributionPlot!(τ_A5_unscaled,"purple")
Plots.savefig(projDir*"/figs/tau_unscaled_Stan.pdf")


######################## Task 3, Priors of expected log reaction #############
# prior for theta was
# theta[j] ~ normal(mu + phi*ISKID[j],tau)

prior_adult = mean(μ_0_unscaled) .+ mean(τ_unscaled) .* randn(N)
prior_kid = mean(μ_0_unscaled) .+ mean(ϕ_unscaled) .+ mean(τ_unscaled) .* randn(N)
makeDistributionPlot(prior_adult,"black")
makeDistributionPlot!(prior_kid,"red")
Plots.savefig(projDir*"/figs/priors_Stan.pdf")

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

Plots.savefig(projDir*"/figs/priors_postOverlay_Stan.pdf")


### Compare to old prior
prior_A5 = mean(μ_A5_unscaled) .+ mean(τ_A5_unscaled) .* randn(N)
makeDistributionPlot(prior_adult,"black",ann=false)
makeDistributionPlot!(prior_kid,"red",ann=false)
makeDistributionPlot!(prior_A5,"green",ann=true)
Plots.savefig(projDir*"/figs/priors_compA5_Stan.pdf")



######################## Task 4, posterior prediction #############
#### a) knowing that it is a child

# Posterior predictive sampling:
idx = Int.(ceil.(rand(N).*N))   # choose random indices which will pick from the simulated posterior
# Sample zlogy according to the model:
zlogy_sim_adult = μ[idx] .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
zlogy_sim_kid = μ[idx] .+ ϕ[idx] .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
# Transform to non-standardised and non-log:
y_sim_adult = exp.(zlogy_sim_adult .* logStd .+ logMean)
y_sim_kid = exp.(zlogy_sim_kid .* logStd .+ logMean)

makeDistributionPlot(y_sim_adult,"black",ann=true)
makeDistributionPlot!(y_sim_kid,"red",ann=true)
Plots.savefig(projDir*"/figs/PP_known_Stan.pdf")

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
zlogy_sim_unknown = μ[idx] .+ ϕ[idx] .* areKids .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
y_sim_unknown = exp.(zlogy_sim_unknown .* logStd .+ logMean)
# Compare in plots:
makeDistributionPlot(y_sim_unknown,"blue",ann=true)
histogram!(y_sim_adult, bins=200, normalize=:pdf, alpha=0.15, linealpha=0.1, color="black")
histogram!(y_sim_kid, bins=200, normalize=:pdf, alpha=0.15, linealpha=0.1, color="red")
Plots.savefig(projDir*"/figs/PP_unknown_Stan.pdf")


##############
#### b-Version 2) not knowing that it is a child, set a fixed fraction
# using fraction 0.5
weight = 0.5
postBeingKid = weight
areKids = postBeingKid .>= rand(N)
idx = Int.(ceil.(rand(N).*N))   # choose a random index which will pick from the simulated posterior
zlogy_sim_unknown = μ[idx] .+ ϕ[idx] .* areKids .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
y_sim_unknown = exp.(zlogy_sim_unknown .* logStd .+ logMean)
makeDistributionPlot(y_sim_unknown,"blue",ann=true)
histogram!(y_sim_adult, bins=200, normalize=:pdf, alpha=0.15, linealpha=0.1, color="black")
histogram!(y_sim_kid, bins=200, normalize=:pdf, alpha=0.15, linealpha=0.1, color="red")
plot!(ann=(1200,0.0025,"Kids: $(weight*100) %"),grid=false)
Plots.savefig(projDir*"/figs/PP_unknown_fixed05_Stan.pdf")



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
    zlogy_sim_unknown = μ[idx] .+ ϕ[idx] .* areKids .+ randn(N).*τ[idx] .+ randn(N).*σ[idx]
    y_sim_unknown = exp.(zlogy_sim_unknown .* logStd .+ logMean)

    makeDistributionPlot!(y_sim_unknown,"blue",ann=false,offset=(i-1)*0.006,scale=1.0)
    makeDistributionPlot!(y_sim_adult,"black",ann=false,offset=(i-1)*0.006,scale=1.0)
    makeDistributionPlot!(y_sim_kid,"red",ann=false,offset=(i-1)*0.006,scale=1.0)
    plot!(ann=(myXlims[2]-200,(i-1)*0.006+0.003,"Kids: $(weight*100) %"),grid=false,yticks=false)
end
# Show plot
plt

Plots.savefig(projDir*"/figs/CompMixture_Stan.pdf")
