Dockerized Delly
================

[Delly](https://github.com/dellytools/delly) is available as a [Docker image](https://hub.docker.com/r/dellytools/delly/). You can pull the Delly image using

`docker pull dellytools/delly`

You can then run Delly from that image. Below we assume your alignment files in BAM format are in /var/data which is mounted as /root in the docker container.

`docker run -it -v /var/data/:/root dellytools/delly`

The container by default starts a shell which you can then use to run Delly.

`delly call -o /root/sv.bcf -g /root/ref.fa /root/control.bam /root/tumor.bam`

Once Delly is finished you can just exit the Docker image.

`exit`
