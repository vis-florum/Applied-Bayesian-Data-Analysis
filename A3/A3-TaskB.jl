using PyCall
using PyPlot
using SpecialFunctions
# uncomment for non-jupyter use:
#pygui(true)

### FUNCTIONS, used for several subtasks:
# %%

# the Bernoulli PDF:
function bernoulli_pdf(θ, y)
    if 0 <= θ <= 1 && (y==0 || y==1)    # Only admissible values
        return θ^y * (1-θ)^(1-y)
    else
        return 0
    end
end

# Likelihood fct of a (series of) coin tosses with bias:
function likelihood_fct(θ,y)
    # Make array s.t. multiple theta values can be input:
    likelihood = Array{Float64}(undef,length(θ))
    # Do that for every theta value requested:
    for ii in 1:length(θ)
        # Elementwise application of the beronoulli_pdf over outcome vector y, then MULTIPLY together:
        likelihood[ii] = prod(bernoulli_pdf.(θ[ii], y))
    end
    return likelihood
end

#%%

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
    likelihood = Array{Float64}(undef,length(θ))
    # Do that for every theta value requested:
    for ii in 1:length(θ)
        # Elementwise application of the log_beronoulli_pdf over outcome vector y, then SUM together:
        likelihood[ii] = sum(log_bernoulli_pdf.(θ[ii],y))
    end
    return likelihood
end

#%%

### Choose and run the subtask:

###########
### Task 2
# Given the parameter, plot outcome probabilities
#%%
θ = 0.5
outcomes = [0, 1]

figure()
bar(outcomes, bernoulli_pdf.(θ,outcomes))
xticks([0, 1])
xlabel("y")
ylabel("θ")
plt.show()

#%%

θ = 0.25

figure()
bar(outcomes, bernoulli_pdf.(θ,outcomes))
xticks([0, 1])
xlabel("y")
ylabel("θ")
plt.show()

#%%

###########
### Task 3
# Given the outcome, plot parameter likelihoods
#%%
θ = 0:0.01:1;
outcome = 1;

figure()
plot(likelihood_fct.(θ,outcome))
plt.show()

#%%
outcome = 0

figure()
plot(likelihood_fct.(θ,outcome))
plt.show()
#%%

###########
### Task 4
# Evaluate likelihood fct after n flips

#%%
θ = 0.5

# a)
println("a) Regular Likelihood:")

for i in 1:3
    n = 10 * 100^(i-1)

    outcome = (rand(n) .>= 0.5) * 1

    print("L(y_1,y_2,...,y_$n | θ = $θ) = ")
    println(likelihood_fct(θ,outcome))
end

println("\nWith rising n, the likelihood -> 0")

#%%
# b) and c)
println("\nb) log-likelihoods, l, and L = exp(l):")

for i in 1:4
    n = 10 * 100^(i-1)

    outcome = (rand(n) .>= 0.5) * 1
    l = log_likelihood_fct(θ,outcome)
    L = exp.(l)

    print("l(y_1,y_2,...,y_$n | θ = $θ) = ")
    println(l)
    print("L(y_1,y_2,...,y_$n | θ = $θ) = ")
    println(L)
end

println("\nYes, n > 100000 could be evaluated without over-/underflow in logarithmic space.")

#%%
# d)
println("\nd) Plotting the likelihood function for a series")
θ = 0:0.01:1;

outcomeSeries = [[1], [1 1], [1 1 0 1]]

for thisOutcome in outcomeSeries
    figure()
    subplot(1,2,1)
    plot(θ,likelihood_fct(θ,thisOutcome),label="L($thisOutcome|θ)")
    grid()
    legend()
    subplot(1,2,2)
    plot(θ,log_likelihood_fct(θ,thisOutcome),label="l($thisOutcome|θ)")
    grid()
    legend()
    plt.show()
end
#%%

###########
### Task 5
# Implement the beta function and create plots

### i)
# Logarithmic version of the beta-distribution
function Log_Beta(θ,a::Any,b::Any)
    #using SpecialFunctions
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


#%%
### ii) sub-tasks 1-3
function Log_posterior(θ,a,b,y)
    heads = sum(y)
    tails = length(y) - heads
    posterior = Array{Float64}(undef,length(θ))
    # Log addition, the brackets are for converting into a single float:
    for ii in 1:length(θ)
        posterior[ii] = (Log_Beta(θ[ii],a,b) .+ log_likelihood_fct(θ[ii], y) .-
                    log(SpecialFunctions.beta(a+heads,b+tails)))[1];
                    #log(factorial(a+heads-1) * factorial(b+tails-1) / factorial(a+heads+b+tails-1)))[1];
    end
    return posterior
end

#%%

println("\nd) Plotting the likelihood function for a series")
θ = 0:0.001:1;

#%%
# Flat prior:
a = 1;
b = 1;

figure()
subplot(1,2,1)
plot(θ,exp.(Log_Beta.(θ,a,b)), label="beta(θ,$a,$b)")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_Beta.(θ,a,b), label="log(beta(θ,$a,$b))")
grid()
legend()
plt.show()

# %%
y = [1]
a = 1
b = 1

figure()
subplot(1,2,1)
plot(θ,exp.(Log_posterior(θ,a,b,y)), label="posterior for $y")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_posterior(θ,a,b,y), label="log-posterior for $y")
grid()
legend()
plt.show()

# Alternatively directly by Beta:
heads = sum(y)
tails = length(y) - heads

figure()
subplot(1,2,1)
plot(θ,exp.(Log_Beta.(θ,a+heads,b+tails)), label="beta(θ,$(a+heads),$(b+tails))")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_Beta.(θ,a+heads,b+tails), label="beta(θ,$(a+heads),$(b+tails))")
grid()
legend()
plt.show()


# %%
y = [1 1]
a = 1
b = 1

figure()
subplot(1,2,1)
plot(θ,exp.(Log_posterior(θ,a,b,y)), label="posterior for $y")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_posterior(θ,a,b,y), label="log-posterior for $y")
grid()
legend()
plt.show()

# Alternatively directly by Beta:
heads = sum(y)
tails = length(y) - heads

figure()
subplot(1,2,1)
plot(θ,exp.(Log_Beta.(θ,a+heads,b+tails)), label="beta(θ,$(a+heads),$(b+tails))")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_Beta.(θ,a+heads,b+tails), label="beta(θ,$(a+heads),$(b+tails))")
grid()
legend()
plt.show()

# %%
y = [1 1 0 1]
a = 1
b = 1

figure()
subplot(1,2,1)
plot(θ,exp.(Log_posterior(θ,a,b,y)), label="posterior for $y")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_posterior(θ,a,b,y), label="log-posterior for $y")
grid()
legend()
plt.show()

# Alternatively directly by Beta:
heads = sum(y)
tails = length(y) - heads

figure()
subplot(1,2,1)
plot(θ,exp.(Log_Beta.(θ,a+heads,b+tails)), label="beta(θ,$(a+heads),$(b+tails))")
grid()
legend()
subplot(1,2,2)
plot(θ,Log_Beta.(θ,a+heads,b+tails), label="beta(θ,$(a+heads),$(b+tails))")
grid()
legend()
plt.show()

#%%

println("The differences are a scaling in the non-log curves and thus
a shifting in the log-curves.
The differences arise because before we only used the LIKELIHOOD function
before, but now we use a normalised POSTERIOR.")

# %%
### iii) More informative prior and reproduce Fig 6.4

# First Column
θ = 0:.01:1
a = 100
b = 100

y = [ones(1,17) zeros(1,3)]

subplot(3,1,1)
plot(θ,exp.(Log_Beta.(θ,a,b)))
grid()
subplot(3,1,2)
plot(θ,exp.(log_likelihood_fct(θ,y)))
grid()
subplot(3,1,3)
plot(θ,exp.(Log_posterior(θ,a,b,y)))
grid()
plt.show()
#%%

# Second Column
θ = 0:.01:1
a = 18.25
b = 6.75

y = [ones(1,17) zeros(1,3)]

subplot(3,1,1)
plot(θ,exp.(Log_Beta.(θ,a,b)))
grid()
subplot(3,1,2)
plot(θ,exp.(log_likelihood_fct(θ,y)))
grid()
subplot(3,1,3)
plot(θ,exp.(Log_posterior(θ,a,b,y)))
grid()
plt.show()

#%%

# Third Column
θ = 0:.001:1
a = 1
b = 1

y = [ones(1,17) zeros(1,3)]

subplot(3,1,1)
plot(θ,exp.(Log_Beta.(θ,a,b)))
grid()
subplot(3,1,2)
plot(θ,exp.(log_likelihood_fct(θ,y)))
grid()
subplot(3,1,3)
plot(θ,exp.(Log_posterior(θ,a,b,y)))
grid()
plt.show()
