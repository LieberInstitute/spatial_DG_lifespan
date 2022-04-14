library("rsconnect")

## Locate app_dir. Edit as needed
app_dir <- here::here("Documents", "code","app_dir")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
load(file.path(app_dir, ".deploy_info.Rdata"), verbose = TRUE)

## Authenticate to shinyapps.io
 # rsconnect::setAccountInfo(
 #     name = deploy_info$name,
 #     token = rsconnect::setAccountInfo(name='libd', token='82689F8C5BB51918D8BF9ED5432FD788', secret='HWc1bsbcJ88HvzGxucEiqcwfGBoTE7ed4x98wC8X'),
 #     secret = deploy_info$secret
 # )
  rsconnect::setAccountInfo(
      name='libd',
      token='82689F8C5BB51918D8BF9ED5432FD788',
      secret='HWc1bsbcJ88HvzGxucEiqcwfGBoTE7ed4x98wC8X'
  )

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
        "spe_wrapper.rds",
        gsub(file.path(app_dir, "www"), "www", dir(file.path(app_dir, "www"), full.names = TRUE))
    ),
    appName = "spatialLIBD_Human_Lymph_Node_10x",
    account = "libd",
    server = "shinyapps.io"
)