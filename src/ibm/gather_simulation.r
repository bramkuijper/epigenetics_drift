# gathers all simulations in certain folder and 
# puts data in one big dataframe


# first some function declarations
# then 

######### first, define functions that do the work #########

# collects parameters from a single file 
# and returns them as a dataframe
#
# Parameters
# ----------
#   filename: str
#       the filename from which parameters 
#       should be collected
#   line_from: int
#       the line from the listing of parameters starts
#   sep:    str
#       separator 
#
# Returns
# -------
#   a dataframe with parameters
collect.params <- function(filename, line_from, sep=";")
{
    # read the table with parameters
    raw.params <- read.table(filename
            ,header= F
            ,sep=sep
            ,skip=line_from # start at the line given by function argument
            ,stringsAsFactors=F) # make sure strings are read as strings

    # what you will get is something like this:
    # V1, V2
    # par1, 0.5
    # par2, 0.9
    # par3, 5

    # we however want
    # par1    par2    par3
    # 0.5     0.9     5
    
    # transpose second column of data frame
    # so that values are lined up in a row 
    p = as.data.frame(t(raw.params$V2), stringsAsFactors=F)

    # now give the rows the names of the columns
    colnames(p) <- raw.params$V1

    # return the data frame
    return(p)
}


# find out where the parameter listing starts
# so that we can read in the data part of the file 
# without having it messed up by the 
# subsequent parameter listing
#
# Parameters
# ----------
#   filename:   str
#       the filename from which the parameters
#       should be read
#
# Returns
# --------
#   linenumber: int/NA
#       returns a number that gives the last line
#       of the data before parameters begin
#       if no last line can be found return NA
last_line_data <- function(filename) 
{
    # read file fil
    f <- readLines(filename)

    # make a reverse sequence
    # of line numbers 
    seqq <- seq(length(f),1,-1)

    # go through all the line numbers
    # of the file (in reverse)
    # and find first line
    # where data is printed (i.e., 
    # a line which starts with a digit)
    for (line_i in seqq)
    {
        # find the first line starting with a digit
        # this is the last line of the data frame
        if (length(grep("^\\d",f[[line_i]])) > 0)
        {
            if (line_i == length(f))
            {
                return(NA)
            }
            else
            {
                return(line_i)
            }
        }
    }

    # no line found? return NA
    return(NA)
}



# function getting all the data from
# a bunch of files in the same directory

# Parameters
# ----------
#   directory: str
#       the directory in which the simulation
#       output files are located
#   sim.file.pattern: regex
#       a pattern of simulation files
#   n.lines.last: int
#       minimum number of lines you want from each simulation
#       (counting from the end that is)
#       if you want the last line from each sim: 1
#       if you want the last 20 lines from each sim: 20
#       if you want all the lines: -1
#
# Returns
# -------
#   dataframe with all the simulation data
gather.sims <- function(directory = "."
             ,sim.file.pattern="sim_.*[[:digit:]]"
             ,n.lines.last=-1)
{
    # get a list of all the simulation files
    all.simulation.files <- list.files(
            path=directory
            ,pattern=sim.file.pattern)

    # place holder variable for a big 
    # data frame with all simulations
    big.dataframe.all.sims <- NULL

    # loop through all the files in the list of files
    for (i in seq(1,length(all.simulation.files)))
    {
        file_i <- all.simulation.files[[i]]

        # filename might be a factor so let's change it to character
        file_i_chr <- as.character(file_i)

        # get the last line of the data of this simulation
        last_line <- last_line_data(filename=file_i_chr)

        # check whether last.line == NA
        # in that case skip
        if (is.na(last_line))
        {
            # cannot find parameters, this simulation
            # did not finish properly and is thus useless
            # continue the next iteration of the loop
            # and forget about this one

            next
        }

        # collect the parameters and put them in a data frame
        param.data <- collect.params(
                                     filename=file_i_chr
                                     ,line_from=last_line + 2)

        # get the actual data and put it in a data frame
        simulation.data <- read.table(
                                      file=file_i
                                      ,header=T
                                      ,sep=";"
                                      ,stringsAsFactors=F
                                      ,nrow=last_line - 1)

        nrow.sim <- nrow(simulation.data)

        data.select.range <- NULL

        if (n.lines.last == -1 || nrow.sim - n.lines.last < 0)
        {
            data.select.range <- seq(1,nrow.sim,1)
        } else {
            data.select.range <- seq(nrow.sim - n.lines.last, nrow.sim,1)
        }

        # modulate simulation data to account for
        # lines one wants to get out of it
        simulation.data <- 
            simulation.data[data.select.range,]

        # merge the parameters and data together
        all.data <- cbind(param.data
                          ,simulation.data
                          ,row.names=NULL)

        all.data[,"file"] <- file_i_chr
        all.data[,"sim.id"] <- i

        big.dataframe.all.sims <- rbind(big.dataframe.all.sims, all.data)
    } # end for

    return(big.dataframe.all.sims)
} # end gather.sims

### here I call the function that I defined before
### it gets the data of all simulations in the current directory
### if you want simulations to be excluded, either put them in a 
### different directory, delete them or select a subset of data
### below
data <- gather.sims(directory=".")

str(data)

