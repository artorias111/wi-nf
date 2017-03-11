# wi-nf

The variant calling pipeline for wild isolates.

### Setup

Information regarding the set of fastqs belonging to a particular isotype are organized in a json document called `isotype_set.json`. This file can be generated using `setup.nf` as follows:

```
nextflow run setup.nf
```

### Usage

```
nextflow run main.nf -resume -e.test=false
```

### Installing Telseq

```
https://gist.githubusercontent.com/danielecook/1ba4db9959cd39641857/raw/fb7bb67952e32e54669e0f64abba7fddc2205708/telseq.rb
```
