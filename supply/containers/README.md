# Description

## Building container for local execution
If you want to run TOGA2 on a local machine or a multiple CPU cluster node, consider building your container with the `apptainer.def` definition file. This defition contains recipes for pulling the latest public TOGA2 version and installing it alongside with all the necessary dependencies. Building a container from this definition file is straightforward and does not require any additional operations:
```
apptainer build <your_container_name> apptainer.def
```
Note that the resulting container supports only local execution. If you want to use TOGA2 with a cluster job manager as Nextflow executor, proceed to the next section.
Running the resulting container is comparably straightforward; all you need is to provide the necessary bindings for directories on the local machine:
```
apptainer exec --bind <your bindings go here>  --env NXF_HOME="<your tmp dir on cluster>/.nextflow" <your_container_name> toga2.py <your_settings>
```
>[!NOTE]
>Nextflow requires a directory ($NFX_HOME) to write temporary files to. Since the suggested container is read-only, the directory must be bound to a local directory on your machine. A good candidate us the user's $TMPDIR:
```
apptainer exec --bind ${TMPDIR}/<other bindings>  --env NXF_HOME="${TMPDIR}/.nextflow" <your_container_name> toga2.py <your_settings>
```

## Building container for cluster execution
Configuring containers for cluster job manager compatibility requires additional setup, including adding directory bindings at build time. The exact configuration, therefore, depends on your batch manager of choice and directory tree structure.

The container recipe provided here, `apptainer_slurm.def`, is a an example build for Slurm compatibility. Note the following important lines:
* `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/host/usr/lib64`: directory `/host/usr/lib64` is a binding alias for local `/usr/lib64` library or a custom local directory containing the necessary libraries. This alias is added to avoid conflict between `/usr/lib64` instances on local machine and within the container.
* The following lines:
    ```
    ## create a dummy user for Slurm
    adduser --disabled-password --gecos "" slurm

    %runscript
        ## by default the container is configured for SLURM; modify the passwd
        echo “slurmadmin:x:300:300::/opt/slurm/slurm:/bin/false” >> /etc/passwd
        echo “slurmadmin:x:300:” >> /etc/group
    ```
    are necessary for Slurm configuration inside the container.
* EXEC_DIR, SLURM_CONF, and LIB_DIR

At runtime, you also have to bind the directories containing Slurm executables and the libraries they use. A sample binding line looks as follows:
```
--bind /opt/slurm,/opt/slurm/bin,/opt/slurm/etc,/var/log/munge,/var/run/munge,/usr/lib64:/host/usr/lib64
```
, as well as export the necessary environment variables:
```
EXEC_DIR=/opt/slurm/bin,SLURM_CONF=/opt/slurm/etc/slurm.conf,LIB_DIR=/host/usr/lib64
```

>[!IMPORTANT]
>Location of the necessary libraries might differ on your machine depending on your file systemt and your batch manager of choice

### Confuring TOGA2 run at runtime
An obvious limitation of running TOGA2 from within the container is code availability: If run from the container, scripts used in parallel steps are inaccessible to your batch manager. One possible way to circumvent this limitation is to run every parallel batch from its own container instance:

<img src="https://github.com/hillerlab/TOGA2/blob/develop/wiki/container_parallel_run.png">

In this case, TOGA2 needs a path to the container from which the parallel jobs will be executed, as well as necessary settings and bindings for parallel container invocation. These can be provided with the following options (added in v2.0.7):
* `--container_image`: a path to the container image to run the code from. This can be any container starting from TOGA2 directory, including containers built from any of the two recipes supplied.
* `--container_executor`: a program to execute the container image; must be available in `$PATH` under the same name. Currently only `apptainer` is available for container execution.
* `--bindigs`: a list of bindings necessary for container invocation. Usually this includes paths to input/output directories and library directories required by your batch manager. The option expects all the bindings as a single comma-separated list, with aliases separated by a colon symbol: `"dir1:alias1,dir2,dir3:alias3"`

>[!IMPORTANT]
>For proper path resolution by both Apptainer and your batch manager, avoid relative paths and symbolic links.

Below is a sample command to run TOGA2 from the container image at `/home/toga2.sif`, with input files stored at `/home/input/` and output being written to `/home/output/`:
```
apptainer exec \
    --bind /home,/opt/slurm,/opt/slurm/bin,/opt/slurm/etc,/var/log/munge,/var/run/munge,/usr/lib64:/host/usr/lib64 \
    --env TMPDIR=/scratch_local,EXEC_DIR=/opt/slurm/bin,SLURM_CONF=/opt/slurm/etc/slurm.conf,LIB_DIR=/host/usr/lib64 \
    /home/toga2.sif toga2.py run \
    --ref_2bit /home/input/ref.2bit \
    --query_2bit /home/input/query.2bit \
    --chain_file /home/input/chains.chain \
    --ref_annotation /home/input/ref.bed \
    --isoform_file /home/input/isoforms.tsv \
    --u12_file /home/input
    -o /home/output \
    -c slurm \
    --container_image /home/toga2.sif \
    --bindings "/home/input,/usr/lib64:/host/usr/lib64"
```