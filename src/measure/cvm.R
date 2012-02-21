library("copula")

N <- c(seq(from = 10, to = 90, by = 10), seq(from = 100, to = 10000, by = 100))
for (n in N) {
    sim <- indepTestSim(n, p = 2, N = 100, print.every = -1)
    content <- paste(sim$dist.global.statistic.independence, collapse = ", ")
    cat("    {", content, "},\n", sep = "")
}
