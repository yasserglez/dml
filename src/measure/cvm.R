library("copula")

cat("{\n")
for (n in seq(from = 100, to = 1000, by = 100)) {
    sim <- indepTestSim(n, p = 2, N = 100, print.every = -1)
    content <- paste(sim$dist.global.statistic.independence, collapse = ", ")
    cat("    {", content, "},\n", sep = "")
}
cat("};\n")
