# neuriteScan
An R Package for neurite image analysis. neuriteScan quantifies neurite density and other significant measures for experiments in Neuroscience area of research.

## Dependencies

neuriteScan uses functions from R Package EBImage and Python modules skimage and numpy via R Package reticulate.

- Install R Package EBImage
```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EBImage")
```

- Install R Package reticulate
```{r}
install.packages("reticulate")
library(reticulate)
py_install("numpy==1.26.4")
py_install("skimage==0.23.1")
py_install("skimage.restoration")
```

## Instalation 

Install neuriteScan from github

```{r}
if (! requireNamespace("devtools", quietly = TRUE))
install.packages("devtools")
devtools::install_github("sysbiolab/neuriteScan")
```
