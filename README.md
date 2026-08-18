<p align="center">
  <p align="center">
    <img width=200 align="center" src="./wiki/LOGO_TOGA2.png" >
  </p>

  <p align="center">
    <img width=200 align="center" src="./wiki/hillerlab.png" >
  </p>

  <span>
    <h1 align="center">
        TOGA2
    </h1>
  </span>

  <p align="center">
    <a href="https://github.com/hillerlab/TOGA2" reference="_blank">
      <img alt="GitHub License" src="https://img.shields.io/github/license/hillerlab/TOGA2?color=blue">
    </a>
  </p>

  <p align="center">
    <samp>
        <span> TOGA2: A faster, more versatile successor of Tool to infer Orthologs from Genome Alignments </span>
        <br>
        <span> The Hiller Lab at the Senckenberg Research Institute </span>
        <br>
        <br>
        <a href="https://github.com/hillerlab/TOGA2/wiki">docs</a> .
        <a href="https://github.com/hillerlab/TOGA2?tab=readme-ov-file#Installation">install</a> .
        <a href ="https://www.biorxiv.org/content/10.64898/2026.06.30.735536v1">preprint</a> .
        <a href="https://hillerlab.com/">us</a> 
    </samp>
  </p>

</p>

---

## Changelog
### v2.0.9
* `run` mode:
    * Paralog annotation is now possible even for transcripts with existing orthologous chains with `--annotate_paralogs` flag. Paralogs annotated this way are retained only if they have loss status of at least *UL*, and no more than 5 paralogs with predicted orthologs per query locus are kept in the annotation.
    * `tmp/` directory is now created at user's `$TMPDIR` addressed to by a local symlink. The original behavior can be restored by setting `--local_tmp` flag.
* `prepare-input` mode:
    * Annotation and isoform file can be now generated from a GFF3/GTF file;
    * Minimal intron length and nonsense-mediated decay (NMD) filters for reference transcripts;
    * Minimal CDS intron length filter;
    * Minimal internal (non-terminal) CDS exon length filter;
    * Reference gene coordinates in BED format (`${ref}.toga.genes.bed`) now added to output;
* `integrate` mode:
    * Adding GTF format query annotation (`query_annotation.gtf.gz`) to output;
    * Adding integration summary (`summary.txt`) to output;
* **NEW MODE**: `orthogroups` enables orthogroup modelling and generates gene family files compatible with [CAFE5](https://github.com/hahnlab/CAFE5) input.
* **NEW MODE**: `toga2agora` creates ordered lists of 1:1 orthologs across multuple runs compatible with the [AGORA](https://github.com/DyogenIBENS/Agora) tool for gene order and synteny analysis.
* **NEW COMPANION PIPELINE**: `toga2kbpython` for estimating gene expression with `kb-python` and TOGA2 output (`v2.0.9d`; early access).
* **NEW COMPANION SCRIPT**: `toga2stats` for assembly quality assessment based on the ancestral gene set.
* Companion dataset:
    * Adding input annotation for two turtle (*Emydura macquarii macquarii*, *Emys orbicularis*) and five percomoprh fish (*Gasterosteus aculeatus*, *Hippoglossus stenolepis*, *Melanotaenia boesemani*, *Nothobranchius furzeri*, *Synchiropus splendidus*) reference species.

For the full list of code changes, see [changelog.md](https://github.com/hillerlab/TOGA2/blob/main/assets/changelog/changelog.md) .

## Installation

### from Github
TOGA2 Makefile provided in this repository contains directives for code compilation, third party software installation, Python package download, and model training. 
```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
cd TOGA2
make
```
>[!IMPORTANT]
> Executing the Makefile <u>will not</u> install Nextflow and Cargo on your machine.

Since the Make directive also installs Python packages globally, you might want to create a dedicated virtual environment beforehand (Python version 3.13 or higher): 
```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
cd TOGA2
python3 -m venv toga2
source toga2/bin/activate
make
```

### from Github + Conda
As an alternative to Python virtual environment, you can also use Conda for environment creation. TOGA2 comes with a Conda configuration file which creates an environment with Python, Cargo, and Nextflow:
```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
cd TOGA2
conda env create -f assets/venv/conda.yaml
conda activate toga2
make
```

### via Apptainer
Container image definition file for Apptainer is provided with TOGA2 under `supply/containers/apptainer.def`. To build the container, make sure you have Apptainer installed, then run the following commands:
```bash
git clone --recurse-submodules https://github.com/hillerlab/TOGA2
TMPDIR=${tmp_dir} apptainer build ${container.sif} supply/apptainer.def
```
where:
* `${tmp_dir}` is the path to the directory to store Nextflow temporary files in;
* `${container.sif}` is the path to save the resulting image to.

The resulting container has `toga2.py` as an entry point. If you run the following or similar command:
```bash
TMPDIR=${tmp_dir} apptainer run --bind ${bound_dir1},${bound_dir2} ${container.sif} toga2.py
```
you should see the TOGA2 start menu.

>[!NOTE]
> The image provided in `supply/` directory contains the latest TOGA2 release, third-party software used for input preparation and TOGA2 annotation, and Nextflow for parallel process management. The container, however, does <ins>not</ins> contain any Nextflow-compatible parallel job executor. To set up your container for batch manager compatibility, see README at `supply/containers` and example recipe at `supply/containters/apptainer_slurm.def`

---

## Running TOGA2
If activated without additional arguments (`toga2.py`), the following start screen is displayed in the user's terminal:
```
Usage: toga2.py [OPTIONS] COMMAND [ARGS]...

  MMP""MM""YMM   .g8""8q.     .g8"""bgd      db          `7MMF'`7MMF'
  P'   MM   `7 .dP'    `YM. .dP'     `M     ;MM:           MM    MM  
       MM     dM'      `MM dM'       `     ,V^MM.          MM    MM  
       MM     MM        MM MM             ,M  `MM          MM    MM  
       MM     MM.      ,MP MM.    `7MMF'  AbmmmqMA         MM    MM  
       MM     `Mb.    ,dP' `Mb.     MM   A'     VML        MM    MM  
     .JMML.     `"bmmd"'     `"bmmmdPY .AMA.   .AMMA.    .JMML..JMML.

  TOGA2 - Tool for Ortholog Inference from Genome Alignment

Options:
  -V, --version      Show the version and exit.
  -help, -h, --help  Show this message and exit.

Commands:
  cookbook            List example commands for 'run' mode
  from-config         Run TOGA2 pipeline with a configuration file
  integrate           Prepare an integrated TOGA2 annotation  by combining annotation with different references
  merge               Merge complementing TOGA2 results for the same reference and query
  postoga             Run postprocessing analysis with Postoga
  prepare-input       Prepare reference annotation files for TOGA2 input
  run                 Run TOGA2 pipeline with command line arguments
  sequence-alignment  Align orthologous sequences from multiple TOGA2 results
  spliceai            Generate SpliceAI predictions for query assembly
  test                Test TOGA2 with companion dataset
```
Except for `toga2.py test`, invoking any of the listed commands without arguments also displays the help message. You can also invoke help message for TOGA2 or any of its commands with `--help/-h` option.

---

## Test run
To ensure that TOGA2 was installed properly, run the following command:

```bash
toga2.py test
```
Provided the successful installation, you will see the TOGA2 execution log printed to stdout, with the results stored at TOGA2/sample_output.
>[!NOTE]
> If you are running TOGA2 from Apptainer, you might have to provide a custom output directory with `--output/-o` option to bypass the read-only container configuration:
```
TMPDIR=${tmp_dir} apptainer run --bind ${bound_dir1},${bound_dir2} ${container.sif} toga2.py test -o ${output_dir}
```

