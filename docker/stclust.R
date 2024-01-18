suppressMessages(suppressWarnings(library('spatialGE')))

# args from command line:
args <- commandArgs(TRUE)


# note, -k --kclusters can be integer or character ('dtc') for valid options
option_list <- list(
    make_option(
        c('-f', '--input_file'),
        help='Path to the count matrix input.'
    ),
    make_option(
        c('-c', '--coordinates_file'),
        help='Path to the barcode spatial coordinates input.'
    ),
    make_option(
        c('-s', '--sample_name'),
        help='Sample name'
    ),
    make_option(
        c('-n', '--normalization'),
        help='Normalization method of `log` or `sct`'
    ),
    make_option(
        c('-k', '--kclusters'),
        help="Number of clusters or `dtc` for auto clusters"
    ),
    make_option(
        c('-o','--output_file_prefix'),
        help='The prefix for the output file'
    )
)

opt <- parse_args(OptionParser(option_list=option_list))

# Check that the file was provided:
# Checks are incomplete
if (is.null(opt$input_file)){
    message('Need to provide a count matrix with the -f/--input_file arg.')
    quit(status=1)
}

if (is.null(opt$coordinates_file)){
    message('Need to provide a count matrix with the -c/--coordinates_file arg.')
    quit(status=1)
}

if (is.null(opt$output_file_prefix)) {
    message('Need to provide the prefix for the output file with the -o arg.')
    quit(status=1)
}

# Set up expected directory structure for spatial data import
indir = paste("./", opt$sample_name, sep="")
dir.create(
    paste(
        "./",
        opt$sample_name,
        "/spatial/",
        sep=""
    ),
    recursive=T
)
# Copy files into directory structure
file.copy(opt$input_file, paste("./", opt$sample_name, sep=""))
file.copy(opt$coordinates_file, paste("./", opt$sample_name, "/spatial/", sep=""))

# Import input data into spatialGE object
spat <- STlist(
    rnacounts=indir, 
    samples=c(opt$sample_name)
)

# Transform the data
spat <- transform_data(spat, method=opt$normalization)

# Run the clustering algorithm
# Note: ws was not included in the parameters
#   it defines the weighting on the distance measure by the clustering algo
#   range of 0-1 expected
#   optional inclusion?
spat <- STclust(
    spat,
    ks=opt$kclusters,
    ws= 0.025
)

# Export of the cluster data
df <- data.frame(
    spat@spatial_meata[opt$sample_name]
)[, c(1,2,3,6)]
colnames(df) <- c("barcodes", "ypos", "xpos", "clusterid")
write.table(
    df,
    paste(opt$output_file_prefix, "tsv", sep="."),
    sep="\t", quote=F, row.names=F
)