library(iCellR)

# Read data into object
LN.data <- load10x("LN/filtered/")

# Make iCellR object
LN.obj <- make.obj(LN.data)

# QC on obj
LN.obj <- qc.stats(LN.obj)
# plot UMIs, genes and percent mito all at once and in one plot. 
# you can make them individually as well, see the arguments ?stats.plot.
stats.plot(LN.obj,
           plot.type = "three.in.one",
           out.name = "UMI-plot",
           interactive = FALSE,
           cell.color = "slategray3", 
           cell.size = 1, 
           cell.transparency = 0.5,
           box.color = "red",
           box.line.col = "green")
