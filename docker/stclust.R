suppressMessages(suppressWarnings(library('spatialGE')))
suppressMessages(suppressWarnings(library('optparse')))

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

# STlist expects that the first column is the gene names- so we don't use row.names arg
rnacounts <- read.table(opt$input_file, sep='\t', header=T, check.names=F)

# Same as for the counts, the expectation is that the coordinates file does not
# have row.names and instead has the barcodes in the first column. For now, however,
# we set the row names and then later alter.
spotcoords <- read.table(opt$coordinates_file, sep='\t', row.names=1, header=T, check.names=T)

# only take the first two columns for the (x,y) positions. Additional columns
# can cause problems downstream
spotcoords <- spotcoords[,c(1,2)]

# the barcodes in coords dataframe can be a superset of the count matrix columns.
# For example, if the matrix is filtered for QC, there may be poor quality spots
# that were filtered out. 
# The opposite is not the case since we cannot have a barcode without a position.
barcodes_from_counts <- colnames(rnacounts)[2:dim(rnacounts)[2]] 
diff_set <- setdiff(barcodes_from_counts, rownames(spotcoords))
num_invalid_barcodes <- length(diff_set)
if (num_invalid_barcodes > 0) {
    max_print <- 5
    if (num_invalid_barcodes < max_print) {
        invalid_barcodes <- paste(diff_set[1:num_invalid_barcodes], collapse=', ')
    } else {
        invalid_barcodes <- sprintf('%s, and %d others.', paste(diff_set[1:max_print], collapse=', '), num_invalid_barcodes-max_print)
    }
    message(sprintf('The set of barcodes in your count matrix must be a subset of those in your coordinates file. Problems include: %s', invalid_barcodes))
    quit(status=1)
} else {
    # the STList constructor below will not accept coordinate files which are a superset of the 
    # count matrix barcodes. We handle that here:
    spotcoords <- spotcoords[barcodes_from_counts,]

    # we need the barcodes in the first col, not the row names. We used the rownames for indexing
    # convenience above, but need to change that now:
    spotcoords <- cbind(rownames(spotcoords), data.frame(spotcoords, row.names=NULL))
}

# Now, to avoid any unexpected issues downstream, we need to conver the column names, etc.
# to preserve the barcodes/column names, we create a dataframe of the original and 'R mutated'
# names. We then run through everything with the mutated names and finally map back.
orig_col_names <- colnames(rnacounts)
proper_names <- make.names(orig_col_names)
colname_mapping = data.frame(
    orig_names = orig_col_names,
    row.names=proper_names,
    stringsAsFactors=F)
colnames(rnacounts) <- proper_names


spotcoords[,1] <- make.names(spotcoords[,1]) 
#rownames(spotcoords) <- make.names(rownames(spotcoords))

# We will use a list of dataframes in the call to STlist
rnacounts_list <- list()
rnacounts_list[[opt$sample_name]] <- rnacounts
spotcoords_list <- list()
spotcoords_list[[opt$sample_name]] <- spotcoords

# Import input data into spatialGE object
spat <- STlist(
    rnacounts=rnacounts_list,
    spotcoords=spotcoords_list, 
    samples=c(opt$sample_name)
)

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