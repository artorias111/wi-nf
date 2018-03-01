#!/usr/bin/bash
# Run the following script to create the virtualenvs required for the pipeline.

pyenv virtualenv 2.7.11 vcf-kit
pyenv local vcf-kit
pip install cython numpy
pip install vcf-kit
pip install https://github.com/AndersenLab/bam-toolbox/archive/0.0.3.tar.gz

pyenv virtualenv 3.6.0 multiqc
pyenv local multiqc
pip install multiqc