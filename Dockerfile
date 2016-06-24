FROM r-base:latest
MAINTAINER Upendra Devisetty <upendra@cyverse.org>
LABEL Description "This Dockerfile is for edgeR multifactorial"

# Run updates
RUN apt-get update
RUN apt-get upgrade -y
RUN apt-get install -y git
RUN apt-get install r-base-dev -y
RUN apt-get install libxml2 -y
RUN apt-get install libxml2-dev -y
RUN apt-get -y install libcurl4-gnutls-dev
RUN apt-get install -y libssl-dev

RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("edgeR");'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("DESeq2");'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("genefilter");'
RUN Rscript -e 'install.packages("devtools", dependencies = TRUE);'
RUN Rscript -e 'devtools::install_github("PF2-pasteur-fr/SARTools")'
RUN Rscript -e 'install.packages("getopt", dependencies = TRUE);'
RUN Rscript -e 'install.packages("knitr", dependencies = TRUE);'
RUN Rscript -e 'install.packages("RColorBrewer", dependencies = TRUE);'

# Add multiple custom functions for Pie_compare, Pie_plot and Bar_compare plots
ENV EDGERM https://github.com/upendrak/edgeR_multi.git
RUN git clone $EDGERM

WORKDIR /edgeR_multi

# change permissions to the wrapper script
ENV BINPATH /usr/bin
RUN chmod +x run_edgeR_multi.r && cp run_edgeR_multi.r $BINPATH
RUN rm -r test_data Dockerfile
RUN mv *.* /

ENTRYPOINT ["run_edgeR_multi.r"]
CMD ["-h"]

# Building and testing
# sudo docker build -t"=ubuntu/iuta:1.0" .
# sudo docker run ubuntu/iuta:1.0 -h
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/iuta:1.0 --gtf mm10_kg_sample_IUTA.gtf --bam1 bam_1 --bam2 bam_2 --fld empirical --test.type SKK,CQ,KY --numsamp 3 --output IUTA_test_1 --groups 1,2 --gene.id Pcmtd1 
# With no gene id (all genes compressed)
# sudo docker run --rm -v $(pwd):/working-dir -w /working-dir ubuntu/iuta:1.0 --gtf mm10_kg_sample_IUTA.gtf --bam1 bam_1 --bam2 bam_2 --fld empirical --test.type SKK,CQ,KY --numsamp 3 --output IUTA_test_1 --groups 1,2
