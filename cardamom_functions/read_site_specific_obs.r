
###
## Default function to read in site specific information for assimilation
###

read_site_specific_obs <- function(variable,infile) {

    # read in the data assuming hearder are present and row do not have numbers
    input=read.csv(infile, header=TRUE)
    # check the variable exists in file
    if (length(which(names(input) == variable)) > 0) {
	# read from the table the desired informatin
	output=input[,variable]
#if (variable == "Evap") {
#    tmp1=input[,"soilevap"]
#    tmp2=input[,"wetevap"]
#output = output - tmp1 - tmp2
#}
    } else {
	# if not then return missing value variable
	output=-9999
    }
    # return the desired variable
    return(output)

} # end of function
