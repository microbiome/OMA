# Create TreeSE object
source("make_TreeSE.R")

# Compare read- and assembly based signals
d <- data.frame(read=colSums(assay(tse, "abundances")),
            assembly=colSums(assay(altExp(tse, "assembly"), "abundances")))


library(ggplot2)
p <- ggplot(d, aes(x=read, y=assembly)) + geom_point(size=1, alpha=0.5) +
       scale_x_log10() +
       scale_y_log10()        
print(p)



