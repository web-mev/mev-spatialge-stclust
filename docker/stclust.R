suppressMessages(suppressWarnings(library('spatialGE')))
suppressMessages(suppressWarnings(library('optparse')))
source('/usr/local/bin/prep_stlist.R')

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
    ),
        make_option(
        c('-x', '--xpos_col'),
        help='The column header for the x-position coordinate metadata'
    ),
    make_option(
        c('-y', '--ypos_col'),
        help='The column header for the y-position coordinate metadata'
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

# transform the name of the normalization scheme:
if (is.null(opt$normalization)){
    message('Need to provide a normalization scheme with the -n/--normalization arg.')
    quit(status=1)
} else if(tolower(opt$normalization) == 'sctransform'){
    norm_scheme <- 'sct'
} else if(tolower(opt$normalization) == 'log'){
    norm_scheme <- 'log'
} else {
    message('We only accept `log` or `SCTransform` for the normalization scheme.')
    quit(status=1)
}

# we can accept an integer or 'automatic' as an argument to STclust below.
# HOWEVER, an integer represented as a string 
if(tolower(opt$kclusters) == 'automatic'){
    cluster_k <- 'dtc'
} else {
    # if here, the argument was not equal to dtc. Can now either be 
    # a number (represented as a string variable) OR an invalid string
    # which cannot be cast as an integer
    as_number <- as.numeric(opt$kclusters)
    if (is.na(as_number)){
        message('The -k/--kclusters option must be an integer or "dtc" for automatic selection.')
        quit(status=1)
    } else {
        cluster_k <- as.integer(as_number)
    }
}

# change the working directory to co-locate with the counts file:
working_dir <- dirname(opt$input_file)
setwd(working_dir)

# call the utility function which will return a list with the necessary items:
spat_list <- prep_stlist(opt$input_file,
                         opt$coordinates_file,
                         opt$sample_name,
                         opt$xpos_col,
                         opt$ypos_col)

# unpack:
colname_mapping <- spat_list$colname_mapping
spat <- spat_list$spat

# Transform the data
spat <- transform_data(spat, method=norm_scheme)

# Run the clustering algorithm
# Note: ws was not included in the parameters
#   it defines the weighting on the distance measure by the clustering algo
#   range of 0-1 expected
#   optional inclusion?
spat <- STclust(
    spat,
    ks=cluster_k,
    ws= 0.025
)

# Export of the cluster data
df <- data.frame(
    spat@spatial_meta[[opt$sample_name]]
)[, c(1,2,3,6)]
colnames(df) <- c("barcodes", "ypos", "xpos", "clusterid")

# convert back to the original barcodes:
df$barcodes <- colname_mapping[df$barcodes, 'orig_names']

if (is.null(opt$output_file_prefix)) {
    output_filename <- sprintf('%s/%s.spatialge_clustered.%s_normed.tsv', working_dir, opt$sample_name, opt$normalization)
} else {
    output_filename <- sprintf('%s/%s.%s.spatialge_clustered.%s_normed.tsv', working_dir, opt$sample_name, opt$output_file_prefix, opt$normalization)
}
write.table(
    df,
    file=output_filename,
    sep="\t", quote=F, row.names=F
)
json_str = paste0('{"clustered_positions":"', output_filename, '"}')
output_json <- paste(working_dir, 'outputs.json', sep='/')
write(json_str, output_json)