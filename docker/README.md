Dockerized DELLY
================

[Delly](https://github.com/dellytools/delly) is available as a [Docker image](https://registry.hub.docker.com/u/dellytools/delly/). You can pull the DELLY image using

`docker pull dellytools/delly`

You can then run DELLY from that image. Below we assume your bam files are in /var/data which is mounted as /data in the docker container:

`docker run -it -v /var/data:/data dellytools/delly /bin/bash`

`delly call -t DEL -o /data/del.bcf -g /data/ref.fa /data/s1.bam /data/s2.bam`

Once DELLY is finished you can just exit the Docker image:

`exit`
