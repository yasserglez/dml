library("CDVine")

family <- 3
basename <- "copula_clayton"
parameters <- c(1e-4, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7,
    8, 9, 10, 20, 30, 40, 50)

values <- seq(from = 0, to = 1, by = 0.1)
u <- rep(values, 11)
v <- as.vector(sapply(values, function (x) rep(x, 11)))

for (parameter in parameters) {
    ## sink(paste(basename, "pdf.dat", sep = "_"), append = TRUE)
    ## r <- BiCopPDF(u, v, family, parameter)
    ## cat(paste(format(r, scientific = TRUE), collapse = "\n"), "\n")
    ## sink(NULL)
	
    ## sink(paste(basename, "cdf.dat", sep = "_"), append = TRUE)
    ## r <- BiCopCDF(u, v, family, parameter)
    ## cat(paste(format(r, scientific = TRUE), collapse = "\n"), "\n")
    ## sink(NULL)

    ## sink(paste(basename, "h.dat", sep = "_"), append = TRUE)
    ## r <- BiCopHfunc(u, v, family, parameter)$hfunc2
    ## r[u == 0] <- 0; r[u == 1] <- 1
    ## cat(paste(format(r, scientific = TRUE), collapse = "\n"), "\n")
    ## sink(NULL)
	
    ## sink(paste(basename, "hinv.dat", sep = "_"), append = TRUE)
    ## r <- .C("Hinv2", as.integer(family), as.integer(121),
    ##        as.double(u), as.double(v),
    ##        as.double(parameter), as.double(0),
    ##        as.double(rep(0, 121)), PACKAGE="CDVine")[[7]]
    ## r[u == 0] <- 0; r[u == 1] <- 1
    ## cat(paste(format(r, scientific = TRUE), collapse = "\n"), "\n")
    ## sink(NULL)
}
