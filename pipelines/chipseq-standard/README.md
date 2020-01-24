# ChIP-Seq Standard Pipeline

User manual for the pipeline can be found here:

https://github.com/NYU-BFX/hic-bench_documentation/blob/master/HiC-Bench_manual.pdf.

ChIP-Seq pipeline usage follows the same structure and format as the HiC pipleline instructions listed there.

This guide assumes that you have cloned this repository to the location `~/hic-bench`. 

## Overview

```bash
# update hic-bench
cd ~/hic-bench
git pull

# create a new project directory
~/hic-bench/code/code.main/pipeline-new-analysis chipseq-standard </path/to/project_directory>

# set input files
cd /path/to/<project_directory>
./code/setup-sample-files.sh <fastq_source_dir>

# create sample sheet template from samples
cd inputs
./code/create-sample-sheet.tcsh <hg19|mm10> <fragment-size>

# fill the control, group, and chip columns in the sample sheet

# run the pipeline
cd /path/to/<project_directory>
./run.dry
./code.main/pipeline-execute <project_name> <name@domain.edu>
```

# Full Walkthrough

## 0. Create a new analysis

If a new analysis project has not already been created, do so with the following command:

```bash
~/hic-bench/code/code.main/pipeline-new-analysis chipseq-standard </path/to/project_directory>
```

Unless otherwise stated, all following steps should be excuted from within the `<project_directory>` you have just created. 

## 1. Set input files

There are two methods to set up your sample files. Pick one of the following:

### Automatic Method (preferred)

From the project directory directory, run:
```bash
./code/setup-sample-files.sh <fastq_source_dir>
```
This will scan the given directory and create symlinks to any found FASTQs in `./inputs/fastq`. It will also clean up file names if they are in the standard Illumina bcl2fastq structure so they can be recognized by `create-sample-sheet.tcsh` (next step). It can be run multiple times to scan multiple directories, but outputs may get overwritten if the sample names are the same. Found FASTQs are printed to keep track of what is happening.

NOTE: Certain sample names could cause issues with this, so be sure to check the output. 

### Manual Method 

Within the corresponding `inputs/fastq` or `inputs/bam` directory, subdirectories should be created with the name of each sample to be included in the analysis. The following naming scheme is preferable:

`<Cell_line>-<ChIP>-<treatment>-<SampleID>`

Each subdirectory should contain all fastq or bam files to be used for that sample through the analysis pipeline. Symlinks can be used if the files are not contained in the same location as the project analysis directory, and are preferable to save storage space. 

## 2. Create project sample sheet

A sample sheet must be created for the analysis project. Run the follow command to do so:

```bash
cd inputs
./code/create-sample-sheet.tcsh <genome> <fragment-size>
```

This scans the `inputs/fastq` directory. If you need to modify the sample names, modify them directly in the `inputs/fastq` first. The sample names in the sample sheet must match the sample names in `inputs/fastq`.

The `<genome>` is `hg19`, `mm10`, etc.

The `fragment-size` should be a numeric argument such as `300`, representing the library size of the sequencing sample. 
Typical fragmentation sizes are 150 for histone marks, and 400 for transcription factors.

A related topic is peak type. This is not specified in the sample sheet, but will be helpful in understanding the results. Both narrow and broad peaks are output by default. Generally, narrow peaks are acceptable for TFs and broad for histone marks. However, exceptions are described in [ENCODE Histone ChIP-seq Data Standards](https://www.encodeproject.org/chip-seq/histone/):

![peak type](https://user-images.githubusercontent.com/6363505/41428326-b6f8e516-6fd7-11e8-9d9d-1c885e01562d.png)

After creation of the sample sheet (should be in `inputs/sample-sheet.tsv`), a manual review process is required to match the correct control or input samples with experimental samples, verify proper grouping names, files, and other entries. 

If editing in Excel, clean up the sample sheet to get rid of problematic characters that get introduced:

```bash
mac2unix sample-sheet.tsv
sed -i 's/"//g' sample-sheet.tsv
```

Return to the main project directory:

```bash
cd ..
```

## 3. Pipeline execution

Test the pipeline with:

```bash
./run.dry
```

Run the pipeline with:

```bash
./code.main/pipeline-execute <project_ID> <name@domain.edu>
```

The pipeline will submit jobs for the various sub-processes. You can monitor the progress with:

```bash
watch qstat
```

### Confirming successful execution

After all the pipeline processes are no longer running, it's important to check for errors. There should be a `*.e*` file (stderr) in the main project directory that displays a brief summary. Usually, only the last few lines are relevant, so you can scroll to the bottom. If the pipeline did not complete successfully, it will specify which step encountered an error.

Each step generates output in `pipeline/<step>`. If there was an error, there is a file written in the `errors` sub-directory. That will tell you which sample had a problem, but the specific error message may not be helpful.

In the `results` sub-directory, there are additional sub-directories for each sample. Within each, the `job.err` file will have the error with relevant context.

## 4. Compile Report

A report template has been supplied in the `report` directory. It will automatically scan the pipeline output, and generate a PDF which includes sample sheets and figures found in the report. 

First, the project info text file needs to be updated with the correct parameters for the project. 

Edit the `project_info.txt` file:
```bash
nano project_info.txt
```

Compile the report:
```bash
cd report
./compile_report_wrapper.sh chipseq_report.Rnw
```

## 5. Copy Results

Most PI's at NYU have a lab directory where sequencing and analysis results are stored. You can copy the analysis results to their lab directory either with an included script, or manually. 

### Results Directory Prep

Analysis results should be copied to the same dated lab subdirectory that contains the source fastq files for the analysis. If source fastq files came from multiple analyses, you should use your discretion to either copy the results to one of them, or create a new lab subdirectory with the date of the analysis. 

So, for example, if our analysis source data came from the `/data/smithlab/2016-06-08/fastq`, ChIP-Seq analysis results should be copied to `/data/smithlab/2016-06-08/ChIP-Seq`. Create this directory with the `mkdir` command, if it does not exist. 

IMPORTANT: Lab directories and files under the path `/ifs/data/sequence/results/` have special permissions automatically applied to them. Therefore, you should not include any arguments in commands such as `cp` or `rsync` which preserve source file permissions (e.g. `cp -a`, `rsync -a`). 

### Automatic Copy

A script has been included in the pipeline to automatically copy the results of many pipeline steps to a PI's results directory.

```bash
./code/code.main/scripts-copy-chipseq-results.sh /path/to/lab_dir/2016-06-08/ChIP-Seq <project_directory>
```
### Manual Copy

If the script does not properly copy all files, you can sync them manually with `rsync`. An easy way to copy is to sync the entire contents of the `results` directory of a given pipeline step, to a destination directory within the analysis results directory. For example, to copy peaks BED files, you could use a command like this.

```bash
# make the destination dir if it doesn't exist
mkdir -p /path/to/lab_dir/2016-06-08/ChIP-Seq/peaks

# copy only 'peaks.bed', ignore all other files
rsync --dry-run -vhrc pipeline/peaks/results/ /path/to/lab_dir/2016-06-08/ChIP-Seq/peaks/ --include="*/" --include="peaks.bed" --exclude="*"
```

## Notes

Errors encountered during pipeline execution can be viewed with:

```bash
./code.main/pipeline-errors
```

Pipeline analysis results can be removed with:

```bash
./code/clean-all
```

