# gathers all simulations in certain folder and puts in big dataframe


# collects parameters from a single file and returns them as a dataframe
collect.params <- function(filename, line_from, sep=";")
{
    raw.params <- read.table(filename
            ,header= F
            ,sep=sep
            ,skip=line_from
            ,stringsAsFactors=F)

    p = as.data.frame(t(raw.params$V2), stringsAsFactors=F)
    colnames(p) <- raw.params$V1
}


# find out where the parameter listing starts
# so that we can read in the data part of the file 
# without having it messed up by the subsequent parameter listing
find_out_param_line <- function(filename) {

    f <- readLines(filename)

    # make a reverse sequence
    seqq <- seq(length(f),1,-1)

    # go through each line in the data file and find first line
    # where data is printed (i.e., a line which starts with a digit)
    for (line_i in seqq)
    {
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            return(line_i + 1)
        }
    }

    return(NA)
}


# get a list of all the simulation files
all.simulation.files <- list.files(
        path="."
        ,pattern="sim_.*")

# print that list to screen
print(all.simulation.files)

# minimum number of lines you want from each simulation
# (counting from the end that is)
# if you want the last line from each sim: 1
# if you want the last 20 lines from each sim: 20
# if you want everything


# place holder variable for a big 
# data frame with all simulations
big.dataframe.all.sims <- NULL

for (file_i in all.simulation.files)
{
    # filename might be a factor so let's change it to character
    file_i_chr <- as.character(file_i)

    param_line <- find_out_param_line(filename=file_i_chr)

    if (!is.na(param_line))
    {
        collect.params(filename=file_i_chr, line_from, sep=";")
    }
}

