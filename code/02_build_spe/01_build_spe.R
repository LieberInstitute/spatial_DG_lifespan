
# cd /dcs04/lieber/marmaypag/lifespanDG_LIBD001/spatial_DG_lifespan/
suppressPackageStartupMessages(library("here"))
# remotes::install_github("drighelli/SpatialExperiment")
# remotes::install_github("LieberInstitute/spatialLIBD")
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

## Create output directories
dir_rdata <- here::here("processed-data", "02_build_spe")
dir.create(dir_rdata, showWarnings = FALSE, recursive = TRUE)

## Define some info for the samples
sample_info <- data.frame(
    sample_id = c(
        "Br8686", 
        "Br2706",
        "Br3942",
        "Br6023",
        "Br8195",
        "Br1412",
        "Br8667",
        "Br5242"
    )
)
sample_info$subject <- sample_info$sample_id
sample_info$sample_path <-
    file.path(
        here::here("processed-data", "01_spaceranger"),
        sample_info$sample_id,
        "outs"
    )
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
## https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/raw-data/sample_info/Visium_HPC_Round1_20220113_Master_ADR.xlsx
## https://github.com/LieberInstitute/spatial_DG_lifespan/blob/main/raw-data/sample_info/Visium_HPC_Round2_20220223_Master_ADR.xlsx
donor_info <- data.frame(
    subject = c("Br8686", "Br2706", "Br3942", "Br6023", "Br8195", "Br1412", "Br8667", "Br5242"),
    age = c(1.05, 17.94, 47.5, 76.38, 0.31, 15.16, 37.3, 73.9),
    sex = c("M", "M", "M", "M", "M", "F", "F", "M"),
    race = c("EA/CAUC", "EA/CAUC", "EA/CAUC", "EA/CAUC", "AA", "EA/CAUC", "EA/CAUC", "EA/CAUC"),
    diagnosis = "Control",
    rin = c(7.1, 8.1, 7.3, 8.4, 7, 7.9, 6.9, 7.4),
    pmi = c(33.5, 33, 27.5, 17.5, 29, 26, 13.5, 21)
)

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id, 
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE
)
Sys.time()

# [1] "2022-03-04 15:19:02 EST"
# 2022-03-04 15:19:05 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
# 2022-03-04 15:22:15 read10xVisiumAnalysis: reading analysis output from SpaceRanger
# 2022-03-04 15:22:28 add10xVisiumAnalysis: adding analysis output from SpaceRanger
# 2022-03-04 15:22:33 rtracklayer::import: reading the reference GTF file
# 2022-03-04 15:24:35 adding gene information to the SPE object
# 2022-03-04 15:24:35 adding information used by spatialLIBD
# [1] "2022-03-04 15:24:50 EST"

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

## Read in cell counts and segmentation results
segmentations_list <-
    lapply(sample_info$sample_id, function(sampleid) {
        file <-
            here(
                "processed-data",
                "01_spaceranger",
                sampleid,
                "outs",
                "spatial",
                "tissue_spot_counts.csv"
            )
        if (!file.exists(file)) {
            return(NULL)
        }
        x <- read.csv(file)
        x$key <- paste0(x$barcode, "_", sampleid)
        return(x)
    })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
    Reduce(function(...) {
        merge(..., all = TRUE)
    }, segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
    segmentations[segmentation_match, -which(
        colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
    )]
colData(spe) <- cbind(colData(spe), segmentation_info)

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 8748
length(no_expr) / nrow(spe) * 100
# [1] 23.90099
spe <- spe[-no_expr, ]


## For visualizing this later with spatialLIBD
spe$overlaps_tissue <-
    factor(ifelse(spe$in_tissue, "in", "out"))

## Save with and without dropping spots outside of the tissue
spe_raw <- spe

saveRDS(spe_raw, file.path(dir_rdata, "spe_raw.rds"))

## Size in Gb
lobstr::obj_size(spe_raw) 
# 1.651702


## Now drop the spots outside the tissue
spe <- spe_raw[, spe_raw$in_tissue]
dim(spe)
# [1] 27853 38287
## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
    message("removing spots without counts for spe")
    spe <- spe[, -which(colSums(counts(spe)) == 0)]
    dim(spe)
}

lobstr::obj_size(spe) 
# 1.534376

saveRDS(spe, file.path(dir_rdata, "spe.rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
