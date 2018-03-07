#!/usr/bin/bash
# Run the following script to create the virtualenvs required for the pipeline.

set -e
export PYENV_ROOT=/usr/local/var/pyenv
if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi
if which pyenv-virtualenv-init > /dev/null; then eval "$(pyenv virtualenv-init -)"; fi
pyenv install 2.7.14
pyenv virtualenv 2.7.14 vcf-kit
pyenv local vcf-kit
pip install cython numpy
pip install vcf-kit
pip install https://github.com/AndersenLab/bam-toolbox/archive/0.0.3.tar.gz

pyenv install 3.6.0
pyenv virtualenv 3.6.0 multiqc
pyenv local multiqc
pip install multiqc