# Copyright
`batchprimer3` is from You et al. (2008) published in [_BMC Genomics_](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-9-253). Source code is reproduced from the supplementary materials for the online version of the article.
It also provides a compiled version of `primer3`.

The skeleton for this def file is taken from [INFRA](https://forgemia.inra.fr/singularity/microservices/apache2). The base image is from [DockerHub](https://hub.docker.com/r/ubuntu/apache2).

# Build and Run Instructions
Building is easy. Ensure `singularity`/`apptainer` is installed, and then the easiest path forward is to build in r/w (sandbox) mode:

```bash
sudo singularity build --sandbox service_apache.sandbox apache.def
```

Then run the server:

```bash
sudo singularity instance start --writable service_apache.sandbox apache2
```

Access the server from a web browser at the following address: [http://localhost:8181/cgi-bin/batchprimer3.cgi](http://localhost:8181/cgi-bin/batchprimer3.cgi).

You can stop the server with:

```bash
sudo singularity instance stop apache2
```

