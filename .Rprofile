## renv
source("renv/activate.R")

## clustermq setup
options(
    clustermq.scheduler = "slurm",
    clustermq.template = "/path/to/file/below" # if using your own template
)

## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded
if (requireNamespace("workflowr", quietly = TRUE)) {
  message("Loading .Rprofile for the current workflowr project")
  library("workflowr")
} else {
  message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
}


## Manage naming conflicts for the benefit of targets
library(targets)
library(conflicted)
conflict_prefer("select", "dplyr", "AnnotationDbi")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("desc", "dplyr")
