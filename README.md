# TurkeyPlot

## Loading packages

Make sure you have ggplot2, magrittr, and Seurat installed.

```{r package loading, warning=FALSE, message=FALSE}
require(ggplot2)
require(magrittr)
require(Seurat)
# can use SeuratData to load example data - you can use your own or get an example from here
require(SeuratData)
# source the TurkeyPlot.R script
source("TurkeyPlot.R")
```

## Making a Turkey plot with your Seurat object

We use the SeuratData pbmc3k but feel free to use your own here! Needs at least two clusters to work.

```{r make a turkey, warning = FALSE, message = FALSE}
pbmc3k <- SeuratData::LoadData("pbmc3k") # will need to install this if not yet
TurkeyPlot(pbmc3k, group.by = "seurat_annotations")

```


