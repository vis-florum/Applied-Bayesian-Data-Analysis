using NamedArrays

taskNr = 3

####### Code used by several tasks #######
function trafoCSV(filename)
    # transforms a given CSV file with 3 columns into a matrix and the into
    # a named array.
    # The first two columns of the CSV file need to be categorical
    # The third column needs to be numerical

    open(filename) do myFile
        lines = readlines(myFile)

        # Header information
        header = split(lines[1],",");
        rowName, colName = header[1:2];

        colNr = length(header)
        rowNr = length(lines) - 1

        # Empty data matrix
        data = Array{String}(undef, rowNr, colNr)

        # Split and put each line into the data matrix:
        for (idx, line) in enumerate(lines[2:end])
            data[idx, 1:colNr] = split(line,",");
        end

        # Get the categorical values:
        hairCols = unique(data[:,1]);
        eyeCols = unique(data[:,2]);
        # Alternatively:
    #        hairCols = Set(data[:,1]);
    #        eyeCols = Set(data[:,2]);

        # Make a named array:
        namedData = NamedArray(zeros(length(eyeCols),length(hairCols)),
                            (eyeCols,hairCols), (colName, rowName))

        # Populate it with absolute numbers:
        for i in 1:length(data[:,1])
            namedData[data[i,2],data[i,1]] = parse(Float64,data[i,3]);
        end

        # Tranfsorm to relative proportions:
        namedData = namedData ./ sum(namedData);

        return namedData
    end
end
###################################

if taskNr == 1
    thisFile = "HairEyeColor.csv";
    println(trafoCSV(thisFile))

elseif taskNr == 2
    thisFile = "HairEyeColor.csv";
    M = trafoCSV(thisFile);

    # a)
    println("\nMarginal of the hair colors:")
    println(sum(M, dims=1))
    # b)
    println("\nMarginal of the eye colors:")
    println(sum(M, dims=2))
    # c )
    println("\nTotal sum:")
    println(sum(M, dims=:))

elseif taskNr == 3
    thisFile = "HairEyeColor.csv";
    M = trafoCSV(thisFile);

    # a)
    print("\nP(Blue eyes, Blond hair) = ")
    println(M["Blue","Blond"])

    # b)
    print("\nP(Brown eyes) = ")
    println(sum(M["Brown",:]))

    # c)
    # Conditional prob -> reallocate on the summed "universe"
    print("\nP(Brown eyes | Red hair) = ")
    println(M["Brown","Red"] / sum(M["Brown",:]))

    # Allow several categories to survive by using a vector as input
    # d) Make intersection of subspaces
    print("\nP((Red OR Blond) AND (Brown OR Blue) = ")
    overlap = sum(M[["Brown","Blue"],["Red","Blond"]]);
    println(overlap)

    # e) Sum the two subspaces (remove the intersection to avoid double counting)
    print("\nP((Red OR Blond) OR (Brown OR Blue) = ")
    combine = sum(M[["Brown","Blue"],:]) + sum(M[:,["Red","Blond"]]) - overlap;
    println(combine)

    # f) Show that P(x,y) != P(x)P(y)
    println("\nP(Blue eyes, Blond hair) != P(Blue eyes) P(Blond hair)")
    print(M["Blue","Blond"])
    print(" != ")
    println(sum(M["Blue",:]) * sum(M[:,"Blond"]))


end
