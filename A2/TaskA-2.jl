using Plots
pyplot()
#using Random
using Distributions
using Optim


taskNr = 3;

####### Code used for several tasks #######
function myNormal(x,m,s)
    # PDF of Gaussian (watch out for elementwise . !)
    y = 1/sqrt(2*pi*s^2) * exp.(-(x.-m).^2 ./ (2*s^2))

    return y
end

# Parameters:
N = 10000;
mu = 3.4;
sigma = sqrt(3);

# Draw samples from built-in normal distribution
r = rand(Normal(mu,sigma), N);
###################################

if taskNr == 1
    # Draw a histogram, normalised
    histogram(r, normalize = :pdf,
        label = string(N, " Samples"))

    # Draw the PDF
    lo = minimum(r);
    hi = maximum(r);
    symRange = maximum([abs(lo)-mu, abs(hi)-mu]);
    x = range(mu-symRange, mu+symRange, length = 100);
    plot!(x, myNormal(x,mu,sigma),
        w = 2,
        label = "PDF")

    savefig("TaskA-2-a.eps")

elseif taskNr == 2
    # Determine granularity and create x and y:
    dx = .1;
    x = collect(-10:dx:20);
    p = myNormal(x,mu,sigma);

    # Riemann sums of
    # Expected value (weigh elementwise by x):
    E_Riemann = sum(p .* x) * dx;
    # Variance (weigh square-distance elementwise by x):
    variance_Riemann = sum(p .* (x.-mu).^2) * dx;

    # Check closeness:
    println(string("Distance of calculated E from true value: ",
        mu - E_Riemann))
    println(string("Distance of calculated E from sample value: ",
        mean(r) - E_Riemann))
    println(string("Distance of calculated Variance from true value ",
        sigma^2 - variance_Riemann))
    println(string("Distance of calculated Variance from sample value ",
        var(r) - variance_Riemann))
    # Note: var(r) automatically scales with n-1 instead of n

elseif taskNr == 3
    # Parameters:
    N = 10000;
    mu = 0;
    sigma = 1;

    # Draw samples from built-in normal distribution
    x = rand(Normal(mu,sigma), N);

    # Make a log normally distributed vector
    y = exp.(x);

    # Draw a histogram, normalised
    histogram(y,
        normalize = :pdf,
        bins = 0:.1:maximum(y),
        w = .1,     # line width of bin edges
        label = string(N, " Samples"))

    # transform the normal pdf to log-normal pdf:
    pLog(x,m,s) = myNormal(log.(x),m,s) ./ x

    # Draw the PDF
    lo = minimum(y);
    hi = maximum(y);
    y_pdf = collect(lo:.1:hi);
    z = pLog(y_pdf,mu,sigma)
    plot!(y_pdf, z,
        w = 1.5,
        label = "PDF")

    savefig("TaskA-2-d-i.eps")

    # Modes using the points in the PDF:
    @time begin
        pdfMax = maximum(z);
        maxIdx = findfirst(z .== pdfMax);
    end

    modeManual = y_pdf[maxIdx]
    println(string("The manual mode is @ y = ", modeManual))

    # Modes using optimisation:
    myObjective(a) = -pLog(a[1],mu,sigma)    # fix the mean and variance
    @time begin
        optim = optimize(myObjective,[0.0, 0.0],NelderMead());
    end

    modeOptim = optim.minimizer[1];
    println(string("The optimised mode is @ y = ", modeOptim))

end
