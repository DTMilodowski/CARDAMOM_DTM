# Function for creating multiple polygons of confidence intervals for use in MHMCMC
# Luke Smallman 13/02/2014

plotconfidence <- function(y) {
#    if (use_parallel) {require(parallel)}
    # Input needs to have dimensions of dim(time,params)
    # If however dimensions are equal to 3 then dim(param,time,cluster)
    # which must be accounted for this this function

    # number of confidence ranges
    nos_confint=10
    # choose colours
    colour_choices=colorRampPalette(c("pink", "red"))
    colour_choices=colour_choices(nos_confint)
    # confidence ranges
    ci_wanted=seq(0.025,0.975,length.out=nos_confint)
    # split between the upper and lower values needed
    upper_ci_wanted=rev(ci_wanted[((length(ci_wanted)*0.5)+1):length(ci_wanted)])
    lower_ci_wanted=ci_wanted[1:(length(ci_wanted)*0.5)]
    ci_wanted=0 ; rm(ci_wanted)

    # Create x value initial
    if (length(dim(y)) == 3){x=1:dim(y)[2]} else { x=1:dim(y)[1]}
    # Just create a blank plot region with axes first. We'll add to this
    #plot(range(x), range(as.vector(y)), type = "n", ann = FALSE)
    for (z in seq(1,nos_confint)) {
        if (length(dim(y)) == 3) {
	   if (z == 1) {
	      # accounting for cluster dimension
	      for (c in seq(1,dim(y)[3])) {
		    if (c == 1) {
			tmp_var=y[sample(1:dim(y)[1],dim(y)[1]),,c]
		    } else {
			tmp_var=tmp_var+y[sample(1:dim(y)[1],dim(y)[1]),,c]
		    }
		}
	    }
	    CI.U <- apply(t(tmp_var),1,quantile,prob=upper_ci_wanted[z],na.rm=TRUE)
	    CI.L <- apply(t(tmp_var),1,quantile,prob=lower_ci_wanted[z],na.rm=TRUE)
        } else {
           # default
           CI.U <- apply(y,1,quantile,prob=upper_ci_wanted[z],na.rm=TRUE)
           CI.L <- apply(y,1,quantile,prob=lower_ci_wanted[z],na.rm=TRUE)
        }
	# Create a 'loop' around the x values. Add values to 'close' the loop
	X.Vec <- c(x, tail(x, 1), rev(x), x[1])
	# Same for y values
	Y.Vec <- c(CI.L, tail(CI.U, 1), rev(CI.U), CI.L[1])
	# Use polygon() to create the enclosed shading area
	# We are 'tracing' around the perimeter as created above
	polygon(X.Vec, Y.Vec, col = colour_choices[z], border = NA)
    } #  confidence interval
    # tidy before we leave
    gc(reset=TRUE, verbose=FALSE)
} # end of function
## Use byte compile
plotconfidence<-cmpfun(plotconfidence)
