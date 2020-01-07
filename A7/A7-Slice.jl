using Distributed
@everywhere using Random
@everywhere using Statistics
@everywhere using LinearAlgebra  # for using eigen
using Plots
using StatsPlots
using SharedArrays
using StatsBase     # for fitting histograms
using CSV, DataFrames   # for saving results as CSV for later


##### Some local setups
# If run from command line externally:
#projDir = dirname(@__FILE__)

# If run from Jupyter/Hydrogen, maybe change to suit you:
projDir= "/lhome/johhub/Desktop/ABDA/A7"

#####
#%% Input Data #################################################################
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
405, 397, 383, 360, 387, 429, 358, 459, 371, 368, 452, 358, 371]

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

# Indicator for each observation i whether it comes from a child or not:
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

# Attempt numbers of each observation i:
@everywhere x = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 1, 2, 3,
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
@everywhere x_alt = Array{Int64,1}[
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
# @everywhere to make them known to all workers
# Progress bar from Jesper:
@everywhere function progressbar(i,N, msg = "")
    if mod(i,round(N/1000))==0
        print("\e[2K") # clear whole line
        print("\e[1G") # move cursor to column 1
        print(msg * " => " * string(round(i/N*100,sigdigits=4)) * "%")
    end
end

@everywhere function mySliceSampler(pdf_log, x_0, w, N = 2000; burnIn = 1000, m=100, verbose=false, msg="")
    # REQUIRES MODULES:
    # Random
    #
    # INPUT:
    # pdf_log --> log-pdf of the target distribution (a function)
    # x_0 --> inital vector of x (dimension Dx1 or 1xD), where D dimension of pdf
    # w --> typical slice size, for a window (dimension Dx1 or 1xD)
    # m --> integer limiting the slice size to m*w
    # N --> number of sample points
    # burnIn --> number of burn-in samples (optional, results will be discarded)
    #
    # OUTPUT:
    # x_s --> sampled x values (NxD), only keep values after full permutations along dimensions
    # pdflog_x_s --> log-pdf of sampled values (Nx1)

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
        if verbose
            progressbar(ii,N+burnIn,msg)
        end

        # Update the new x progressively for each dimension:
        L   = 1*x_0  # R and L need to be of same dimension, because we apply the pdf upon those!
        R   = 1*x_0
        x_1 = 1*x_0  # "1*" needed, for correct type
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
        x_s[i,:] = x_1     # different directions in arrays n'importe, sample point, after one round of permutations
        pdflog_x_s[i] = pdflog_x_1 # now contains log_pdf of all new coordinates of new point, since has been updated at each permutation
        # note that x_0 has been updated during the permutations
    end
    return x_s, pdflog_x_s
end


@everywhere function mySliceSampler_Cov(pdf_log, x_0, C, N; burnIn = 1000, m=100, verbose=true, msg="")
    # REQUIRES MODULES:
    # Random
    # LinearAlgebra
    #
    # INPUT:
    # pdf_log --> log-pdf of the target distribution (a function)
    # x_0 --> inital vector (mode of a previous run) of x (dimension 1xD), where D dimension of pdf
    # C --> Covariance matrix of a previous burn-in of the same log-pdf
    # m --> integer limiting the slice size to m*w
    # N --> number of sample points
    # burnIn --> number of burn-in samples (optional, results will be discarded)
    #
    # OUTPUT:
    # x_s --> sampled x values (NxD), only keep values after full permutations along dimensions
    # pdflog_x_s --> log-pdf of sampled values

    λ,Eigvec = LinearAlgebra.eigen(C)        # eigenvalue and matrix of eigenvectors, C is the covariance matrix of the coordinates of a previous run
    w = sqrt.(abs.(λ))    # the single standard-deviations along the dimensions

    D = length(x_0)  # the dimension of the distribution (nr of coordinates)
    x_s = Array{Float64}(undef,N,D)   # Samples (NxD), preallocate array, with undefined values
    pdflog_x_s = Array{Float64}(undef,N)   # (unnormalised) log-prob of sampled points (Nx1)
                                           # for reducing evaluations

    pdflog_x_1 = pdf_log(x_0);       # first value is already known!
    wd = zeros(D)       # preallocating window

    L = similar(x_0)
    R = similar(x_0)
    x_1 = similar(x_0)
    # Using dots for array assignment is a bit faster


    for ii in 1:(N+burnIn)
        if verbose
            progressbar(ii,N+burnIn,msg)
        end

        # Update the new x progressively for each dimension:
        L   .= x_0  # R and L need to be of same dimension, because we apply the pdf upon those!
        R   .= x_0
        x_1 .= x_0  # "1*" needed, for correct type
        for d in randperm(D)
            ### 1) Make the Slice
            # random "vertical" position
            #vertical = rand() * pdf_log(x_0) # in normal space
            z = pdflog_x_1 + log(rand())                # saving evaluations, pdflog_x_1 was updated for current permutation in the end

            ### 2) Make the window, "Stepping-out" alogrithm
            # Randomly place a window, containing x_0:
            # Here the entire boundaries are updates simultaneously, to
            # account for the covariation between parameters!
            wd = w[d]*Eigvec[:,d]    # associated co-std value times its eigenvector gives estimate of window
            L .= x_0 .- rand() .* wd
            R .= L .+ wd

            # Randomly share the max window size among left/right:
            J = floor(m*rand());
            K = (m-1) - J;

            # Extend window to the left, until outside slice or allowance seizes:
            while ((J > 0) && (z < pdf_log(L)))
                #println("Lefting")
                L .-= wd
                J -= 1
            end

            # Extend window to the right, until outside slice or allowance seizes:
            while ((K > 0) && (z < pdf_log(R)))
                #println("Righting")
                R .+= wd
                K -= 1
            end

            ### 3) Sample from window
            # finding an allowable point:
            while true
                # Updating over all dimensions simultaneously, to account for the
                # covariation between parameters:
                x_1 .= L .+ rand() .* (R .- L)

                # this + breaking out of loop reduces the amount of log-pdf evaluations:
                pdflog_x_1 = pdf_log(x_1)

                if (pdflog_x_1 >= z)     # new value found
                    x_0 .= 1*x_1     # update value of all dimensions
                    break
                end
                # Value was not within slice, shrink the interval:
                #elseif all(x_1 .> x_0)
                if sign.(R .- x_0) == sign.(x_1 .- x_0)
                    R .= 1*x_1
                else
                    L .= 1*x_1
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

# Diagonal Slice sampler using burn-ins, adapted from Jesper's code:
@everywhere function diagonalSliceSampling(pdf_log, x_0, w, N = 10_000, burnIn = 1_000; m = 100, verbose = true)
    # pre-burn-in run for getting a reasonable start point for the burn-in:
    x_s, l_pdf = mySliceSampler(pdf_log, x_0, w, burnIn; burnIn=0, m=m, verbose=verbose, msg="pre-burn-in")
    x_0 = x_s[argmax(l_pdf),:]   # Mode of pre-burn-in is new start point
    w = std(x_s, dims=1)  # a new guess for the window, my samples are column vectors

    println()
    # burn-in run for getting an idea about parameter covariation (bread-crums):
    x_s, l_pdf = mySliceSampler(pdf_log, x_0, w, burnIn; burnIn=0, m=m, verbose=verbose, msg="burn-in")

    println()
    # actual diagonal slice Sampling run, taking advantage of parameter Covariation
    C = cov(x_s)
    x_0 = x_s[argmax(l_pdf),:] # Mode of burn-in is new start point
    return mySliceSampler_Cov(pdf_log, x_0, C, N; burnIn=0, m=m, verbose=verbose, msg="MCMC run")
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
# accomodate, that the means no longer are a constant, but an array, i.e.
# one mean for each theta
@everywhere function log_prior_theta(θ::Array{Float64},μ::Array{Float64},σ::Float64)
    # Log-normal distribution
    # Logarithmise the normal pdfs in series:
    # We can leave out the 1/sqrt(2π) since it does not change sampling from the posterior
    return sum(-log(σ) .- 0.5 * ((θ .- μ) ./ σ).^2)
end

@everywhere function log_lklhd_fct(J,jointParams::Array{Float64},zlogy::Array{Float64,1},zx::Array{Float64,1},ind::Array{Int64,1},K_j::Array{Int64,1})
    # return the joint lklhd of the given input data
    # Reparametrise, s.t. we sample eta ~ N(0,1) instead of theta. Theta is then reconstructed.
    η_0 = jointParams[1:J]
    η_1 = jointParams[(J+1):(2*J)]
    μ_0, ϕ_0, μ_1, ϕ_1, σ, τ_0, τ_1 = jointParams[(2*J+1):end]

    # Reconstruct thetas (span up by kids indicator K_j):
    θ_0 = (μ_0 .+ ϕ_0 * K_j) .+ τ_0 * η_0
    θ_1 = (μ_1 .+ ϕ_1 * K_j) .+ τ_1 * η_1

    jointLklhd = 0.0

    # the case of sigma <= 0 is already caught by the posterior main function

    # Logarithmise the normal pdfs in series:
    # summing up:
    # note, that we go through each observation, vectorising is hard:
    # zlogy ~ N(theta0 + theta1*zx, sigma)
    for i in 1:length(zlogy)
        j = ind[i]
        regMean = θ_0[j] + θ_1[j] * zx[i]   # regression mean
        jointLklhd += -log(σ) - 0.5 * ((zlogy[i] - regMean) / σ)^2
    end

    return jointLklhd
end


# Skipping the prior funciton speeds up the process.
# Still needs concretising by data
# Avoid using global vars, since slow
@everywhere function log_posterior_fct(J,jointParams::Array{Float64}, zlogy, zx, ind, child_j)
    # The joint log-pdf posterior of the model
    η_0 = jointParams[1:J]
    η_1 = jointParams[(J+1):(2*J)]
    μ_0, ϕ_0, μ_1, ϕ_1, σ, τ_0, τ_1 = jointParams[(2*J+1):end]

    #μ = μ_0 .+ ϕ .* K_j.*1.0    # current sampled μ for each individual -> this evaluation from A6 has now moved to the likelihood function due to reparametrisation!

    if (σ > 0) && (τ_0 > 0) && (τ_1 > 0)
        return log_lklhd_fct(J,jointParams::Array{Float64}, zlogy, zx, ind, child_j) +
               log_prior_theta(η_0,zeros(J),1.0) +  # eta0 ~ N(0,1)
               log_prior_theta(η_1,zeros(J),1.0) +  # eta1 ~ N(0,1)
               #log_prior_flat(μ_0) +   # is zero
               #log_prior_flat(μ_1) +   # is zero
               #log_prior_flat(ϕ_0) +   # is zero
               #log_prior_flat(ϕ_1) +   # is zero
               log_prior_flat_positiveOnly(σ) +
               log_prior_flat_positiveOnly(τ_0) +
               log_prior_flat_positiveOnly(τ_1)
    else
        return -Inf
    end
end


#%% Logarithmise and Transform Data ############################################
@everywhere J = maximum(ind)    # number of individuals
@everywhere I = length(y)       # number of observations

# Go to log-space:
@everywhere logy = log.(y)

@everywhere logMean = mean(logy)
@everywhere logStd = std(logy)

@everywhere trainMean = mean(x)
@everywhere trainStd = std(x)

# Mean-centre and scale each element:
@everywhere zlogy = (logy .- logMean) ./ logStd;

# Standardise training data:
@everywhere zx = (x .- trainMean) ./ trainStd;

#% Concretise the posterior by including the data #############################
@everywhere log_posterior(jointParams::Array{Float64}) =
            log_posterior_fct(J,jointParams::Array{Float64}, zlogy, zx, ind, child_j)


################################################################################
############################# MCMC #############################################
#%% Chain Setup
noParams = 2*J+7
noOfChains = 4
N = 4*10^6                # total number of samples
N_chain = convert(Int64, N / noOfChains)   # samples per chain
burnIn = 10^4   # the larger the burn-In the better the performance (better diagonal)
w = ones(Float64,noParams) * 0.1         # typical window size
m = 100                             # multiplier for maximum window size

#%% MCMC run, make several chains
chain = SharedArray{Float64}(N_chain,noParams,noOfChains);
chainPdf = SharedArray{Float64}(N_chain,noOfChains);

# Prerun to initialise the functions:
x_start = rand(Float64,noParams) * 0.5
@time trash1, trash2 = diagonalSliceSampling(log_posterior,x_start,w,100,1000,m=m,verbose=true);

# (using threads for parallel computing (other stuff didn't work)):
# Threads.@threads led to problems too
# now using @distributed, prefix with @sync to wait for completion
@sync @distributed for k in 1:noOfChains
    x_start = rand(Float64,noParams) * 0.5   # initial vector
    chain[:,:,k], chainPdf[:,k] = diagonalSliceSampling(log_posterior,x_start,w,N_chain,burnIn,m=m,verbose=true);
end

# My machine: 4 cores, i7-7500U 2,7 GHz x 2 each, on Linux
# N=10^5 takes around
# N=10^6 takes around 9 min
# N=4*10^6 takes around 40 min


################################################################################
############################# RESULTS ##########################################

##############
# Acess results and Transform back (undo mean-centering and scaling and go to non-log space)
# See derivation in PDF file instead

# Access the parameters and fuse the chains:
η_0 = Array{Float64,2}(undef,N,J)
η_1 = Array{Float64,2}(undef,N,J)
θ_0 = Array{Float64,2}(undef,N,J)
θ_1 = Array{Float64,2}(undef,N,J)
for j in 1:J
    η_0[:,j] = chain[:,j,:][:]
    η_1[:,j] = chain[:,J+j,:][:]
end

μ_0 = chain[:,2*J+1,:][:]
ϕ_0 = chain[:,2*J+2,:][:]
μ_1 = chain[:,2*J+3,:][:]
ϕ_1 = chain[:,2*J+4,:][:]
σ   = chain[:,2*J+5,:][:]
τ_0 = chain[:,2*J+6,:][:]
τ_1 = chain[:,2*J+7,:][:]

for j = 1:J
    θ_0[:,j] = (μ_0 .+ ϕ_0 .* child_j[j]) .+ τ_0 .* η_0[j]
    θ_1[:,j] = (μ_1 .+ ϕ_1 .* child_j[j]) .+ τ_1 .* η_1[j]
end


### Exporting
CSV.write(projDir*"/sliceSamples/export_samples_many.csv",
          DataFrame([θ_0 θ_1 μ_0 ϕ_0 μ_1 ϕ_1 σ τ_0 τ_1]),
          writeheader=true)

### Importing, if rerun
# # read into a DataFrame:
# chain_csv = CSV.read(projDir*"/sliceSamples/export_samples.csv"; comment="#", normalizenames=true)
#
# θ_0 = Array{Float64,2}(undef,N,J)
# θ_1 = Array{Float64,2}(undef,N,J)
# for j in 1:J
#     θ_0[:,j] = chain_csv[:,j]
#     θ_1[:,j] = chain_csv[:,J+j]
# end
# μ_0 = chain_csv[:,2*J+1]
# ϕ_0 = chain_csv[:,2*J+2]
# μ_1 = chain_csv[:,2*J+3]
# ϕ_1 = chain_csv[:,2*J+4]
# σ   = chain_csv[:,2*J+5]
# τ_0 = chain_csv[:,2*J+6]
# τ_1 = chain_csv[:,2*J+7]
###


#####
# Autocorrelation and ESS
ess_θ_0 = Array{Float64,1}(undef,J)
ess_θ_1 = Array{Float64,1}(undef,J)
for j in 1:J
    ess_θ_0[j] = ess(θ_0[:,j])
    ess_θ_1[j] = ess(θ_1[:,j])
    println("ESS of θ_0_$j = ",ess_θ_0[j])
    println("ESS of θ_1_$j = ",ess_θ_1[j])
end
println("Minimal ESS of θ_0 is for individual $(argmin(ess_θ_0)) := ",ess_θ_0[argmin(ess_θ_0)])
println("Minimal ESS of θ_1 is for individual $(argmin(ess_θ_1)) := ",ess_θ_1[argmin(ess_θ_1)])
println("Maximal ESS of θ_0 is for individual $(argmax(ess_θ_0)) := ",ess_θ_0[argmax(ess_θ_0)])
println("Maximal ESS of θ_1 is for individual $(argmax(ess_θ_1)) := ",ess_θ_1[argmax(ess_θ_1)])

println("ESS of μ_0 = ",ess(μ_0))
println("ESS of ϕ_0 = ",ess(ϕ_0))
println("ESS of μ_1 = ",ess(μ_1))
println("ESS of ϕ_1 = ",ess(ϕ_1))
println("ESS of σ   = ",ess(σ))
println("ESS of τ_0 = ",ess(τ_0))
println("ESS of τ_1 = ",ess(τ_1))


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
    unlogExpected = Array{Float64,1}(undef,N)
    unlogExpected = exp.(θ_0_unscaled[:,j] .+
                θ_1_unscaled[:,j] .* x .+
                0.5 .* σ_unscaled.^2);
    return unlogExpected
end

function E_y_group(x,isKid)
    # isKid is either 1 or 0
    # transform to zx:
    #zx = (x - trainMean) / trainStd
    unlogExpected = Array{Float64,1}(undef,N)
    unlogExpected = exp.(μ_0_unscaled + ϕ_0_unscaled .* isKid .+
                    (μ_1_unscaled .+ ϕ_1_unscaled .* isKid) .* x .+
                    0.5 .* (σ_unscaled.^2 .+ τ_0_unscaled.^2 .-
                            τ_tilde_1_unscaled.^2 .* x .+ τ_1_unscaled.^2 .* x^2))
    return unlogExpected
end


function curveSwarm(nrCurves,j,N,α,subpl=0)
    checkpoints = 100
    attempts = range(0,21,length=checkpoints)
    curves = Array{Float64,2}(undef,N,checkpoints)

    for (i, att) in enumerate(attempts)
        curves[:,i] = E_y_ind(att,j)
    end

    if subpl > 0
        plot!(legend=false,ylims=[100,800],subplot=subpl)
    else
        plot(legend=false,ylims=[100,800])
    end
    randIdces = Int.(ceil.(rand(nrCurves).*N))

    for i = 1:nrCurves
        if subpl > 0
            plot!(attempts,curves[randIdces[i],:],
                  color="blue",linealpha=α,subplot=subpl,
                  legend=false,ylims=[100,800])
        else
            plot!(attempts,curves[randIdces[i],:],
                  color="blue",linealpha=α,
                  legend=false,ylims=[100,800])
        end
    end
end

function curveSwarmGroup(nrCurves,isKid,N,α,subpl=0)
    checkpoints = 100
    attempts = range(0,21,length=checkpoints)
    curves = Array{Float64,2}(undef,N,checkpoints)

    for (i, att) in enumerate(attempts)
        curves[:,i] = E_y_group(att,isKid)
    end

    #if subpl > 0
    #    plot(legend=false,ylims=[100,800],subplot=subpl)
    #else
    #    plot(legend=false,ylims=[100,800])
    #end
    randIdces = Int.(ceil.(rand(nrCurves).*N))

    for i = 1:nrCurves
        if subpl > 0
            plot!(attempts,curves[randIdces[i],:],
                  color=RGB(isKid,0,0),linealpha=α,subplot=subpl,
                  legend=false,ylims=[100,800])
        else
            plot!(attempts,curves[randIdces[i],:],
                  color=RGB(isKid,0,0),linealpha=α,
                  legend=false,ylims=[100,800])
        end
    end
end


################################################################################
############################# TASKS ############################################

##### Task 3
for j in [1,3,4]
    ### Diagnostics to compare to Jesper:
    makeDistributionPlot(exp.(θ_1_unscaled[:,j]),"blue")
    plot!(xlabel="exp(theta_1[$j])")
    Plots.savefig(projDir*"/figs/th1$j-Slice.pdf")

    ### Expected reaction times:
    makeDistributionPlot(E_y_ind(1,j),"blue")
    makeDistributionPlot!(E_y_ind(5,j),"red")
    plot!(xlims=[200,800],grid=false,xlabel="E[y]-ind$j")
    Plots.savefig(projDir*"/figs/times-$j-Slice.pdf")

    ### Credible regression lines:
    curveSwarm(200,j,N,.15)
    col = RGB(child_j[j],0,0)
    scatter!(1:1:sum(ind.==j),y[ind.==j],
        markersize=10,markercolor=col,markerstrokecolor=col)
    plot!(grid=false,xlabel="attempt nr",ylabel="reaction time")
    Plots.savefig(projDir*"/figs/swarm-$j-Slice.pdf")
end


### All Regression lines:
plot(layout=(4, 9),size=(1500, 1000),
     legend=false,grid=false,xaxis=false,yaxis=false)
for j in 1:J
    curveSwarm(200,j,N,.15,j)
    if (j-1)%9==0
        myax=true
    else
        myax=false
    end
    col = RGB(child_j[j],0,0)
    scatter!(1:1:sum(ind.==j),y[ind.==j],
        markersize=4,markercolor=col,markerstrokecolor=col,
        subplot=j, ann=(12,700,"j = $j"),xaxis=true,yaxis=myax,
        framestyle=:axes)
end
Plots.savefig(projDir*"/figs/allcurves-Slice.pdf")


##### Diagnostics to compare to Jesper:
makeDistributionPlot(τ_0,"black")
Plots.savefig(projDir*"/figs/tau0-Slice.pdf")
makeDistributionPlot(τ_0_unscaled,"black")
Plots.savefig(projDir*"/figs/tau0-unsc-Slice.pdf")

makeDistributionPlot(σ,"black")
Plots.savefig(projDir*"/figs/sigma-Slice.pdf")
makeDistributionPlot(σ_unscaled,"black")
Plots.savefig(projDir*"/figs/sigma-unsc-Slice.pdf")


##### Extra
makeDistributionPlot(exp.(ϕ_1_unscaled .* 1),"black")

##### Swarm of the groups
attempts = range(0,22,length=200)
curve_mean_adult = mean.(E_y_group.(attempts,0))
curve_mean_kid = mean.(E_y_group.(attempts,1))

plot(attempts,curve_mean_adult,
    color=:black,ylims=[100,800],linewidth=8,linealpha=0.6)
plot!(attempts,curve_mean_kid,
      color=:red,ylims=[100,800],linewidth=8,linealpha=0.6)
curveSwarmGroup(200,0,N,.15)
curveSwarmGroup(200,1,N,.15)
plot!(grid=false,xlabel="attempt nr",ylabel="reaction time")
Plots.savefig(projDir*"/figs/swarm-groups-Slice.pdf")


### Chain diagnostics
using LaTeXStrings
chstart = 650000
chend = 655000
plot(chstart:chend,θ_0[chstart:chend,:],legend=false,xlabel="sample nr",ylabel=L"\theta_{0_j}")
Plots.savefig(projDir*"/figs/stuckChain-theta0-Slice.pdf")

plot(chstart:chend,θ_1[chstart:chend,:],legend=false,xlabel="sample nr",ylabel=L"\theta_{1_j}")
Plots.savefig(projDir*"/figs/stuckChain-theta1-Slice.pdf")

plot(chstart:chend,σ[chstart:chend],legend=:topleft,xlabel="sample nr",label=L"\sigma")
plot!(chstart:chend,ϕ_0[chstart:chend],xlabel="sample nr",label=L"\varphi_0")
plot!(chstart:chend,ϕ_1[chstart:chend],xlabel="sample nr",label=L"\varphi_1")
plot!(chstart:chend,μ_0[chstart:chend],xlabel="sample nr",label=L"\mu_0")
plot!(chstart:chend,μ_1[chstart:chend],xlabel="sample nr",label=L"\mu_1")
Plots.savefig(projDir*"/figs/stuckChain-restPars-Slice.pdf")
plot(chstart:chend,τ_0[chstart:chend],xlabel="sample nr",label=L"\tau_0")
plot!(chstart:chend,τ_1[chstart:chend],xlabel="sample nr",label=L"\tau_1")
Plots.savefig(projDir*"/figs/stuckChain-tau-Slice.pdf")

makeDistributionPlot(τ_1,"red")
plot!(grid=false,xlabel=L"\tau_1")
Plots.savefig(projDir*"/figs/tau1-Slice.pdf")
