function getPost(testResult)

    # Create a Dictionary containing the possible combinations:
    P = Dict()

    # Idea: solve by recursive calls of function,
    # constantly updating the prior --> not implemented

    deseaseFrequency = 0.001
    truePositives    = 0.99
    falsePositives   = 0.05

    ### Priors:
    P[":("] = deseaseFrequency
    P[":)"] = 1.0 - P[":("]

    ### Likelihoods
    # Positive test restuls, given status:
    P["+|:("] = truePositives
    P["+|:)"] = falsePositives
    # Negative test restuls, given status:
    P["-|:("] = 1.0 - truePositives
    P["-|:)"] = 1.0 - falsePositives

    # Update for each new measurement:
    for i in 1:length(testResult)
        # Marginal Likelihoods:
        P["+"] = P["+|:("] * P[":("] + P["+|:)"] * P[":)"]
        P["-"] = P["-|:("] * P[":("] + P["-|:)"] * P[":)"]

        ### Posteriors:
        P[":(|+"] = P["+|:("] * P[":("] / P["+"]
        P[":(|-"] = P["-|:("] * P[":("] / P["-"]
        P[":)|+"] = 1.0 - P[":(|+"]
        P[":)|-"] = 1.0 - P[":(|-"]

        # Update priors to the posteriors:
        P[":("] = P[string(":(|",testResult[i])]
        P[":)"] = 1.0 - P[":("]
    end

    # Desired Posteriors are the new Priors:
    # If no test result, then this also works, because the posterior = prior
    posteriors = [P[":("], P[":)"]]

    return posteriors
end

# Enter the received results:
myResult = "+-"

# Show the posteriors for each step when receiving test results,
# including "no result" at the start:
for i in 1:length(myResult)+1
    print("\nAccumulated Results: ")
    println(myResult[1:i-1])
    println("[P(:(|T), P(:)|T)] = ",getPost(myResult[1:i-1]))
    print("Sanity check: Sum(P) = ")
    println(sum(getPost(myResult[1:i-1])))
end
