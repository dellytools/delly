Dockerized DELLY
================

This folder is used to create a Docker image of DELLY. Installation and usage instruction can be found on the [docker](https://www.docker.io) website.
Once you have docker installed you can just pull the DELLY image:

`sudo docker pull trausch/delly`

You can then run DELLY from that image. Below we assume your bam files are in /var/data which is mounted as /data in the docker container:

`sudo docker run -i -t -v /var/data:/data trausch/delly /bin/bash`

`/delly/src/delly -t DEL -o /data/del.vcf -g /data/ref.fa /data/s1.bam /data/s2.bam`

Once DELLY is finished you can just exit the Docker image:

`exit`
