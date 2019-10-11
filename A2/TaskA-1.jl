using Plots
pyplot()

taskNr = 2;

function coinFlips(n,p)
    # Generate a random vector from the uniform distribution
    r = rand(n);
    heads = zeros(Float64,n);

    # Turney operator, make head or tail (1,0)
    heads = map(x -> x >= p ? 1 : 0, r)

    # Calculate the cumulative ratio of heads along the series
    headsNr = cumsum(heads)
    flipnr = collect(1:n)
    headsRatio = headsNr ./ flipnr

    return heads, headsRatio
end

if taskNr == 1
    N = 1000;
    P = 0.5;    # ratio of tails
    (h, hRatio) = coinFlips(N,P);

    # Plot it
    plot(hRatio,
         xlabel = "Flips",
         xscale = :log10,
         xlims  = (1,N),
         ylabel = "Frequency Heads",
         yticks = 0:.1:1,
         marker = :d,
         label  = "")
    # for holding the plot
    plot!([1,N],[P,P],     # start @ 1 to make log work
          linestyle = :dash,
          label = "",
          w = 2)
    annotate!(N/10, abs(hRatio[1] - (P-0.2)), # located either above or below the 0.5 line
             text(string("Final ratio = ", hRatio[end])))

    savefig("TaskA-1-a.eps")
    savefig("TaskA-1-a.png")

elseif taskNr == 2
    for i = 1:3
        N = 100 * 10^(i-1);
        P = 0.75;   # ratio of tails

        (h, hRatio) = coinFlips(N,P);
        histogram(h,
            label = string(N, " flips"),
            bar_width = 1,
            normalize = :probability)
        plot!([0, 0],[0, P],
            label = "PMF",
            linestyle = :solid,
            c = :red,
            w = 3)
        plot!([1, 1],[0, 1-P],
            label = "",
            linestyle = :solid,
            c = :red,
            w = 3)
        scatter!([0, 1],[P, 1-P],
            label = "",
            marker = (0,20,:red))

        savefig(string("TaskA-1-c-", i, ".eps"))
        savefig(string("TaskA-1-c-", i, ".png"))
    end
end
