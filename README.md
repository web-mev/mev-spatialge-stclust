# SpatialGE clustering

This repository contains a WebMeV-compatible tool for running clustering on spatial transcriptomics using the spatialGE R package (https://github.com/FridleyLab/spatialGE). 

To run this tool, pull the Docker container (see https://github.com/web-mev/mev-spatialge-stclust/pkgs/container/mev-spatialge-stclust). Change to the directory where your data lives, then run:

```
docker run -d -v $PWD:/work ghcr.io/web-mev/mev-spatialge-stclust:<tag> \
    Rscript /usr/local/bin/stclust.R \
        -f /work/<COUNTS FILE> \
        -c /work/<COORDS FILE> \
        -s <SAMPLE_NAME> \
        -n <METHOD> \
        -k <NUM_CLUSTERS> \
        -o <OUTPUT_PREFIX>
```
If successful, this should create an output file mapping the barcodes of the spaces to the positions and cluster identifiers..

Note:
- `<tag>`: This is the fully-qualified reference to the Docker image, which is tagged to correspond to the Git commit hash. This ensures each Docker image is clearly linked with a specific state of the repository.
- `/work/<COUNTS FILE>`: This is a path to your tab-delimited count file, with genes in rows and barcodes/spots in columns. Since we assume you have changed to the directory containins your raw data and mounted it in the Docker container at `/work`, we write this as `/work/<COUNTS FILE>`.
- `/work/<COORDS FILE>`: This is a path to your tab-delimited coordinates file which has three columns (in order):
    - The barcodes. These barcodes must match those in the count matrix (the column headers there).
    - The x/horizontal position of the spot.
    - The y/vertical position of the spot. 
As with the counts file, we assume you have changed to the directory containins your raw data and mounted it in the Docker container at `/work`, we write this as `/work/<COORDS FILE>`.
- `<SAMPLE_NAME>` is the name of the sample. We incorporate this name into the name/path of the output normalized expression.
- `<METHOD>` is either `log` or `SCTransform`, verbatim. Case does not matter.
- `<NUM_CLUSTERS>` is either an integer (>1) or "dtc" which enables automatic selection of this parameter. 
- `<OUTPUT_PREFIX>` is an optional prefix which will be added to the front of the output file name.