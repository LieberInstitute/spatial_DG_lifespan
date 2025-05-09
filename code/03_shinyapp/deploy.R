library("rsconnect")
library("here")
library("withr")

## Locate app_dir. Edit as needed
app_dir <- here::here("code", "03_shinyapp")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
#load(file.path(app_dir, ".deploy_info.Rdata"), verbose = TRUE)

source(file.path(app_dir, "token.R"))

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
## Note that appFiles has to be relative to app_dir.
## Drop the www directory if you didn't customize the documentation files and
## edit app.R accordingly.
rsconnect::deployApp(
    appDir = app_dir,
    appFiles = c(
        "app.R",
		"pseudobulk_spe.rds",
		"modeling_results.rds",
		"sig_genes_subset.rds",
		with_dir(here("code", "03_shinyapp"), dir("spe_shiny", full.names = TRUE)),
		with_dir(here("code", "03_shinyapp"), dir("www", full.names = TRUE))
    ),
    appName = "Lifespan_DG",
    account = "libd",
    server = "shinyapps.io"
)
