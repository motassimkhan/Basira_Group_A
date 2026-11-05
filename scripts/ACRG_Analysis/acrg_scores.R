
.libPaths(c(file.path(path.expand("~"), "R_libs"), .libPaths()))
library(estimate); library(MCPcounter)
args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]; outdir <- args[2]

# ESTIMATE
filt <- file.path(outdir, "ESTIMATE_Filtered.gct")
res  <- file.path(outdir, "ESTIMATE_Results.gct")
filterCommonGenes(input.f = inp, output.f = filt, id = "GeneSymbol")
estimateScore(input.ds = filt, output.ds = res, platform = "affymetrix")
sc  <- read.delim(res, skip = 2, row.names = 1, check.names = FALSE)
ids <- setdiff(colnames(sc), "Description")
write.csv(data.frame(GSM = ids,
                     ESTIMATE_Immune_Score  = as.numeric(sc["ImmuneScore", ids]),
                     ESTIMATE_Stromal_Score = as.numeric(sc["StromalScore", ids])),
          file.path(outdir, "ACRG_ESTIMATE.csv"), row.names = FALSE)

# MCP-counter, all populations
expr <- as.matrix(read.delim(inp, row.names = 1, check.names = FALSE))
mcp  <- MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
write.csv(data.frame(GSM = colnames(mcp), t(mcp), check.names = FALSE),
          file.path(outdir, "ACRG_MCPcounter_all.csv"), row.names = FALSE)
cat("done\n")
