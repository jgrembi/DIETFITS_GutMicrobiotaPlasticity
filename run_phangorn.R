# Mape a tree from multiple alignment using phangorn package

# { time Rscript --default-packages=methods --no-save --no-restore --verbose run_phangorn.R > dietPhangorn.Rout 2> errorDietPhangorn.Rout ; } 2> timePhangorn.txt
library(ape)
library(phangorn)
path.to.dir <- "path/to/output"
alignment <- readRDS(file.path(path.to.dir, "decipher_align.rds"))

alignment <- as(alignment, "matrix")
phang.align <- phyDat(alignment, type = "DNA")
dm <- dist.ml(phang.align)
# Note, treNJ tip order != sequence order
treeNJ <- NJ(dm) 
fit <- pml(treeNJ, data = phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic") #, control = pml.control(trace = 0))

fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, 
                    optGamma = TRUE, rearrangement = "stochastic")
detach("package:phangorn", unload=TRUE)

save(list = c("phang.align", "dm", "treeNJ", "fit", "fitGTR"),
     file = file.path(path.to.dir, "phangorn_fit.rds"))


