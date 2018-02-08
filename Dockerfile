FROM linuxbrew/linuxbrew:1.4.3

USER root 
# Install latest version of R

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ed \
        less \
        locales \
        vim-tiny \
        wget \
        ca-certificates \
        fonts-texgyre \
        libxml2-dev \
        libcairo2-dev \
        libsqlite-dev \
        libmariadbd-dev \
        libmariadb-client-lgpl-dev \
        libpq-dev \
        libssh2-1-dev \
    && rm -rf /var/lib/apt/lists/*

## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

## Use Debian unstable via pinning -- new style via APT::Default-Release
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
    && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default

ENV R_BASE_VERSION 3.4.3

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
RUN apt-get update \
    && apt-get install -t unstable -y --no-install-recommends --allow-unauthenticated \
        littler \
                r-cran-littler \
        r-base=${R_BASE_VERSION}* \
        r-base-dev=${R_BASE_VERSION}* \
        r-recommended=${R_BASE_VERSION}* \
        && echo 'options(repos = c(CRAN = "https://cran.rstudio.com/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
    && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
    && ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
    && ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
    && ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
    && install.r docopt \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

USER linuxbrew
# Install homebrew files
RUN brew install gcc \
    && brew install https://raw.githubusercontent.com/Linuxbrew/homebrew-core/043fb1f50af078db481b971d36c605f0dcf72ccd/Formula/jdk.rb \
    && brew tap brewsci/science \
    && brew install \
            bwa \
            samtools \
            bcftools \
            bedtools \
            nextflow \
            sambamba \
            vcflib \
            vcftools \
            python2 \
            picard-tools \
            pigz \
            parallel \
            snpeff \
            muscle \
            vcfanno \
            igvtools \
            bamtools \
            trimmomatic


RUN brew install fastqc --ignore-dependencies

RUN pip2 install numpy cython multiqc
RUN pip2 install https://github.com/AndersenLab/bam-toolbox/archive/0.0.3.tar.gz vcf-kit \
    && ln /home/linuxbrew/.linuxbrew/bin/python2 /home/linuxbrew/.linuxbrew/bin/python

# Take over the R lib
RUN sudo chown -R linuxbrew:linuxbrew /usr/local/
ENV R_LIBS_USER=/usr/local/lib/R/site-library
# Install R packages and link python
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile \
    && Rscript -e 'install.packages(c("memoise", "tidyverse", "cowplot", "ggmap", "ape", "devtools", "knitr", "rmarkdown", "aws.s3"))' \
    && Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("phyloseq"))' \
    && Rscript -e 'devtools::install_github("andersenlab/cegwas")'

# Install telseq
RUN wget -O /home/linuxbrew/.linuxbrew/bin/telseq https://github.com/zd1/telseq/raw/master/bin/ubuntu/telseq  \
    && chmod +x /home/linuxbrew/.linuxbrew/bin/telseq

