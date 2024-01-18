# placeholder container name
FROM ghcr.io/web-mev/base-spatialge-docker

# Install the optparse package, which allows better commandline arg parsing
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN mkdir -p /opt/software
ADD stclust.R /opt/software/