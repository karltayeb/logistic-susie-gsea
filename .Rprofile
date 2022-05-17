## Activate the renv package to manage the project's package library
# source("renv/activate.R")

## clustermq
options(
  clustermq.scheduler = "slurm",
  clustermq.template = "/path/to/file/below" # if using your own template
)

## Attach the libraries we will always need to work in the console
if(interactive())
  suppressPackageStartupMessages(
    {
      library(targets)
      library(workflowr)
    }
  )

# Check the project build status when opening the project
if(interactive())
  message("Run 'wflow_build(\"analysis/m_00_status.Rmd\")' to see the project status")
