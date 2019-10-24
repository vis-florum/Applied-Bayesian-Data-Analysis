### PARALLELISATION:
# start julia with "julia -p auto" or "julia -p 4" to use multiple processors
using Revise
using Random
using SpecialFunctions
using PyCall
using PyPlot
using Plots
using StatsPlots
using CmdStan

#%%

# If run from command line externally:
#projDir = dirname(@__FILE__)

# If run from Jupyter/Hydrogen, maybe change to suit you:
projDir= "/home/johhub/Desktop/ABDA/A4"
#projDir= "/lhome/johhub/Desktop/ABDA/A4"
tmpDir = projDir*"/tmp"

################################################################

function runStanMC(y,a,b,burnIn,N,projDir,tmpDir,keepchains=false,noOfChains=4)
    # REQUIRES MODULES:
    # CmdStan
    #
    # INPUT:
    # y --> observed data
    # burnIn --> number of bur in samples
    # N --> total number of sample points
    # projDir --> directory where the analysis should be conducted, and plots should be saved
    # tmpDir --> directory where the temp files can be stored
    # keepchains (optional) --> keep all vectors of the various chains, false
    #                           per defaul, i.e. all chain are thrown in one vector
    #
    # OUTPUT:
    # x_s --> sampled x values for all chains (Nx4 array)

    ### 1) Define the Stan model, note that variable names must
    # correspond to observed data later
    # dichtonomous model
    dichtModel = "
    data {
    int<lower=0> J; // number of flips in the observed data, setting the minimum amount of data under lower
    int<lower=0,upper=1> y[J];  // coin flips
    }
    parameters {
    real<lower=0,upper=1> theta; // prob of getting a head
    }
    transformed parameters { // no transformed variables to use
    }
    model {
    theta ~ beta($a, $b);         // prior distribution for theta
    y ~ bernoulli(theta);       // likelihood, note that stan will create the posterior automatically.
    }
    ";

    ### 2) Data for the Stan model, note that variable names must
    # correspond to defined model in stan-code
    observedData = Dict("J" => length(y), "y" => y);

    ### 3) Chain specs
    myModel = Stanmodel(
                    Sample(save_warmup=true,
                           num_warmup=burnIn,
                           num_samples=N,
                           thin=1),   # thin: Period between saved samples
                    name = "flipping",
                    model = dichtModel,
                    printsummary = false,
                    tmpdir = tmpDir,
                    nchains = noOfChains);    # number of chains (4 default)

    ### 4) Run CmdStan:
    # rc:     return code (0 if all is fine)
    # chn:    chain results
    # cnames: vector of variable names
    rc, chn, cnames = stan(myModel,
                           observedData,
                           projDir,
                           diagnostics = false,
                           CmdStanDir = CMDSTAN_HOME);

    ### Results without burn-in
    #chns = chn[(n_burnin+1):end, :, :]
    # the 7th (out of 8) coulmn in the .value method are the actual sampled points
    # the 3rd dimension are the different chains: 1,2,3 --> chain 1,2,3

    x_s = 1* chn.value[(burnIn+1):end,7,:]  # the 1* converts the "Axis Array" in the chain into an Array{Float64,2}

    if keepchains == false
        x_s = reshape(x_s,:);   # Throw all the chain samples in a single vetor
    end

    return x_s
end

#%%

# Chain inputs:
N = 10^6
N_burnin = 1000

#%%

# First Column of Figure 6.4:
a = 100 # prior specs
b = 100 # prior specs
y = [ones(Int64,17); zeros(Int64,3)];    # make an Array{Intt64,1} array

# to receive the prior curve, just set the data to empty array []:
chain = runStanMC([],a,b,N_burnin,N,projDir,tmpDir);
histogram(chain, bins=100, normalize=:pdf, label="prior")
density!(chain,linewidth=3,label="prior")

# to receive the likelihood curve, just set the prior to a=1, b=1 (flat):
chain = runStanMC(y,1,1,N_burnin,N,projDir,tmpDir);
histogram!(chain, bins=100, normalize=:pdf, label="likelihood")
density!(chain,linewidth=3,label="likelihood")

posterior1 = runStanMC(y,a,b,N_burnin,N,projDir,tmpDir);
histogram!(posterior1, bins=100, normalize=:pdf, label="posterior")
density!(posterior1,linewidth=3,label="posterior")

Plots.savefig("A4-A-C1_Stan.pdf")

#%%
# Second Column
a = 18.25
b = 6.75
y = [ones(Int64,17); zeros(Int64,3)];    # make an Array{Intt64,1} array

# to receive the prior curve, just set the data to empty array []:
chain = runStanMC([],a,b,N_burnin,N,projDir,tmpDir);
histogram(chain, bins=100, normalize=:pdf, label="prior")
density!(chain,linewidth=3,label="prior")

# to receive the likelihood curve, just set the prior to a=1, b=1 (flat):
chain = runStanMC(y,1,1,N_burnin,N,projDir,tmpDir);
histogram!(chain, bins=100, normalize=:pdf, label="likelihood")
density!(chain,linewidth=3,label="likelihood")

posterior2 = runStanMC(y,a,b,N_burnin,N,projDir,tmpDir);
histogram!(posterior2, bins=100, normalize=:pdf, label="posterior")
density!(posterior2,linewidth=3,label="posterior")

Plots.savefig("A4-A-C2_Stan.pdf")

#%%
# Third Column
a = 1
b = 1
y = [ones(Int64,17); zeros(Int64,3)];    # make an Array{Intt64,1} array

# to receive the prior curve, just set the data to empty array []:
chain = runStanMC([],a,b,N_burnin,N,projDir,tmpDir);
histogram(chain, bins=100, normalize=:pdf, label="prior")
density!(chain,linewidth=3,label="prior")

# to receive the likelihood curve, just set the prior to a=1, b=1 (flat):
chain = runStanMC(y,1,1,N_burnin,N,projDir,tmpDir);
histogram!(chain, bins=100, normalize=:pdf, label="likelihood")
density!(chain,linewidth=3,label="likelihood")

posterior3 = runStanMC(y,a,b,N_burnin,N,projDir,tmpDir);
histogram!(posterior3, bins=100, normalize=:pdf, label="posterior")
density!(posterior3,linewidth=3,label="posterior")

Plots.savefig("A4-A-C3_Stan.pdf")

#%%

density(posterior1, label="a=100, b=100")
density!(posterior2, label="a=18.25, b=6.75")
density!(posterior3, label="a=1, b=1")
Plots.savefig("A4-A-all_Stan.pdf")
