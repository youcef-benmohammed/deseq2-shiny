# DESeq2 Shiny app
#   docker build -t deseq2-shiny .
#   docker run --rm -p 3838:3838 deseq2-shiny   ->  http://localhost:3838
FROM bioconductor/bioconductor_docker:RELEASE_3_19

RUN Rscript -e 'BiocManager::install(c("DESeq2", "shiny", "shinyjs", "shinythemes", "DT", \
      "plotly", "ggplot2", "heatmaply", "pheatmap", "rmarkdown"), ask = FALSE, update = FALSE)'

WORKDIR /app
COPY . /app

EXPOSE 3838
CMD ["Rscript", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
