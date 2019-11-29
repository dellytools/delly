You can build a [delly](https://github.com/dellytools/delly) singularity container (SIF file) using

`sudo singularity build delly.sif delly.def`

Once you have built the container you can run analysis using

`singularity exec delly.sif delly --help`
