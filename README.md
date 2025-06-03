# geneva-hrd
Graaphical interface (R/shiny) to the oncoscanR package for computing the HRD score used by the Geneva HRD test.

https://github.com/yannchristinat/oncoscanR

https://www.bioconductor.org/packages/release/bioc/html/oncoscanR.html

## Installation
First install required R packages:
```
install.packages(c("shiny", "shinycssloaders", "magrittr", "jsonlite", "readr"))
```

Then install oncoscanR package from Bioconductor:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("oncoscanR")
```

## Running the app
Open the file `app.R` in Rstudio and click on "Run App" button.

Alternatively, deploy the app to your local shiny server or shinyapps.io (https://shiny.posit.co/r/articles/share/shinyapps/).

