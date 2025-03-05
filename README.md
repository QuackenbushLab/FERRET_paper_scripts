# Replicating FERRET Paper Results

## Installation and Setup

You will need to install and set up FERRET and the GRN inference methods evaluated in the paper. These instructions assume a Unix/Mac environment with `sudo` access.

### FERRET Installation and Setup

To install and set up FERRET, follow the instructions in the [GitHub repository](https://github.com/QuackenbushLab/FERRET).

### Python Installation and Setup

You will need Python to run scGeneRai. To install Python, run the following commands:
```
sudo apt-get update
sudo apt-get install python3
```
You will also need to install `conda`. To do this, run the following commands:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
export PATH=$PATH:”/home/ubuntu/miniconda3/bin”
```

### MATLAB Installation and Setup

You will need MATLAB to run SINGE and GRISLI. To install MATLAB, do the following steps:
1.  Download MATLAB. You may have institutional access. Institutional access for Harvard is available [here](https://www.mathworks.com/academia/tah-portal/harvard-university-30596681.html).
2.  Download the [license](https://www.mathworks.com/licensecenter), then select Install and Activate.
3.  Run ```sudo unzip -X -K matlab_R2023b_glnxa64.zip -d matlab_2023b_installer```, substituting version number as needed.

For GRISLI, you will also need the Bioinformatics Toolbox and the Statistics and Machine Learning Toolbox. To install these, run:
```
wget https://www.mathworks.com/mpm/glnxa64/mpm
chmod 777 mpm
sudo ./mpm install --release R2024a --products Bioinformatics_Toolbox Statistics_and_Machine_Learning_Toolbox
```

### Julia Installation and Setup

You will need Julia to run PIDC. To install Julia, run the following commands:
```
wget https://julialang-s3.julialang.org/bin/linux/x64/1.10/julia-1.10.0-linux-x86_64.tar.gz
tar zxvf julia-1.10.0-linux-x86_64.tar.gz
vi ~/.bashrc
export PATH="$PATH:/home/ubuntu/julia-1.10.0"
sudo snap install julia --classic
```

### R Installation and Setup

You will need R to run scSGL and LEAP. To install R, run the following commands:

```
sudo apt update
sudo apt install r-base r-base-dev -y
```

### Docker Installation and Setup

You will need to install Docker to run SINGE. To install Docker, run the following commands:
```
sudo apt install docker.io
sudo usermod -aG docker ubuntu
docker login -u <username>
```

### GRISLI Installation and Setup

To install GRISLI, run the following commands:
```
sudo apt-get install libgtk2.0-0
git clone https://github.com/PCAubin/GRISLI
wget https://files.inria.fr/spams/spams-matlab-precompiled-v2.6-2017-02-27-Linux.tar.gz
tar -zxvf spams-matlab-precompiled-v2.6-2017-02-27-Linux.tar.gz
cp -r spams-matlab-v2.6/ GRISLI/
```

### LEAP Installation and Setup

LEAP is available from the [CRAN archive](https://cran.r-project.org/src/contrib/Archive/LEAP/).

### PIDC Installation and Setup

To install PIDC, run `julia` on the command line and issue the following lines of code interactively:
```
> import Pkg
> Pkg.add(“Pkg”)
> Pkg.add("InformationMeasures")
> Pkg.add(“PyPlot”)
> Pkg.add(“LightGraphs”)
> Pkg.add(“GraphPlot”)
> Pkg.add(“NetworkInference”)
```

### ScGeneRai Installation and Setup

1.  To install scGeneRai, run the following commands:

    ```
    git clone https://github.com/PhGK/scGeneRAI.git
    sudo apt install python3-pip
    pip3 install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu
    pip install tqdm 
    cd ../scGeneRAI/
    pip install -r requirements.txt
    pip install pandas
    ```

2.  Add the path to scGeneRai to the PYTHONPATH in `.bashrc` and to the `PATH` in `.profile`.


### scSGL Installation and Setup

1. Run ```git clone https://github.com/Single-Cell-Graph-Learning/scSGL.git```.
2. Run ```cd scSGL```.
3. Check the file requirements in `requirements.txt` and install each one manually using `pip`.
4. Launch `R` in the command line and run ```install.packages("pcaPP")```.

### SINGE Installation and Setup

To install and run SINGE using the Docker container, run the following commands:

```
git clone https://github.com/gitter-lab/SINGE.git
cd SINGE
wget https://hastie.su.domains/glmnet_matlab/glmnet_matlab_new.zip
unzip glmnet_matlab_new.zip -d SINGE/
sudo apt-get update
docker pull agitter/singe:0.5.1
```

### SCORPION Installation and Setup

SCORPION is available from [CRAN](https://cran.r-project.org/web/packages/SCORPION/index.html) and [GitHub](https://github.com/kuijjerlab/SCORPION)

### SCENIC Installation and Setup

SCENIC can be installed from GitHub using the instructions [here](https://rdrr.io/github/aertslab/SCENIC/).

## Preprocessing

In most cases, you should not need to preprocess the data again. Instead, simply download them using the FERRET R package function  ```DownloadFERRETData()```. However, if you wish to replicate the preprocessing steps, the instructions for doing this are below.

### Preprocessing CPTAC Data

1. Download the CPTAC data from [Genomic Data Commons](https://portal.gdc.cancer.gov/analysis_page?app=Projects) using the following steps:
    1. Set the Experimental Strategy to "scRNA-seq".
    2. Select "CPTAC-3" from the list of projects.
    3. Select Biospecimen > TSV and save to your machine.
    4. Select Clinical > TSV and save to your machine.
    5. Select Manifest and save to your machine.
    6. Download the scRNA-seq data using the [GDC Data Transfer Tool](https://docs.gdc.cancer.gov/Encyclopedia/pages/Manifest_File/).
2. Run **CPTAC_singlecell_preprocessing.R**, setting the following variables:
   -  **pythonpath:** The path to your Python installation. This will be needed by the Reticulate R package.
   -  **dataDir:** The path to the location where you downloaded the CPTAC data.
   -  **biospecimenDir:** The name of the local folder (inside dataDir) where you installed the biospecimens
   -  **clinicalDir:** The name of the local folder (inside dataDir) where you installed the clinical information
3. Run **assign_pseudotime_single_cell_cptac.R** to assign pseudotime based on gene expression, setting the following variables:
   -  **expression_dir:** The path to the directory where the gene expression data are stored
   -  **pseudotime_result_dir:** The path to the directory where you wish to store the pseudotime assignment.
4. Run **assign_pseudotime_single_cell_cptac.R** to assign pseudotime based on gene expression, setting the following variables:
   -  **expression_dir:** The path to the directory where the gene expression data are stored
   -  **pseudotime_result_dir:** The path to the directory where you wish to store the pseudotime assignment.
5. Run **assign_pseudotime_single_cell_manual_cptac.R** to manually assign a pseudotime starting point and endpoint based on the gene expression data and a Monocle trajectory, setting the following variables:
   -  **expression_dir:** The path to the directory where the gene expression data are stored
   -  **pseudotime_result_dir:** The path to the directory where you wish to store the pseudotime assignment.
6. Run **single_cell_PC_expression_CPTAC_subset.R** to project CPTAC gene data onto the first 100 principal components, subset the number of cells, and split into folds, setting the following variables:
   -  **pca_result_dir:** The file path where you wish to store the results of PCA
   -  **plot_result_dir:** The file path where you wish to store the PCA plots
   -  **expression_dir:** The path to the file containing the preprocessed expression data
   -  **split_dir:** The path to the directory where you wish to store the data after splitting into folds
7. Run **single_cell_CPTAC_splits_subset_pseudotime_only.R** to run the corresponding splits for the assigned pseudotime trajectory, setting the following variables:
   -  **split_result_dir:** The path to the directory where you wish to store the data after splitting into folds
   -  **expression_dir:** The path to the file containing the preprocessed expression data
   -  **split_result_dir_pseudotime:** The path to the directory where you wish to store the pseudotime data after splitting into folds

### Preprocessing HTAN Data

1. Download the HTAN SCLC data from the [HTAN website](https://humantumoratlas.org/explore) and accompanying metadata using the following filters:
   -  Organ = Lung NOS
   -  Atlas = MSK
   -  Assay = scRNA-seq
   -  File Format = csv
2. Remove the following samples, which do not correspond to lung tissue:
   -  RU 681 and RU 255 (brain metastasis)
   -  RU 666 (bone metastasis)
   -  MRM35429319 and 35379369 (spinal fluid)
   -  1559_1262C P96 and RU1141 (not listed in metadata file)
   -  RU1128 (liver metastasis)
3. Run **preprocess_htan_msk.R**, setting the following variables:
   -  **dataDir:** The path to the directory where the data were downloaded
4. Run *assign_pseudotime_single_cell_htan_msk.R*, setting the following variables:
   -  **dataDir:** The path to the directory where the data were downloaded
5. Run *assign_pseudotime_manual_single_cell_htan_msk.R*, setting the following variables:
   -  **pseudotime_result_dir:** The path to the directory where the initial pseudotime assignments were stored
6. Run *single_cell_HTAN_splits_subset_all_cells.R*, setting the following variables:
   -  **expression_dir:** The path to the file containing the preprocessed expression data
   -  **split_result_dir:** The path to the directory where you wish to store the data after splitting into folds
   -  **pseudotime_dir:** The path to the directory with the final pseudotime assignments
   -  **pseudotime_result_dir:** The path to the directory where you wish to store the pseudotime data after splitting into folds

## ChIP-seq Data Evaluation

### Download Data

1. Navigate to [CISTROME DB](http://dc2.cistrome.org/#/)
2. Select the following filters:
   -  Species = Homo sapiens
   -  Biological Sources = B Lymphocyte, B-lympocytes, fibrobasts (sic), Fibroblast
3. Download only those files meeting the following criteria:
   -  Passing all quality controls
   -  Genome is known and is not hg18
   -  Autosomes are available
   -  The factor is available for both fibroblasts and B-lymphocytes
   -  The cells are untreated
4.  Download the GENCODE genomes for [hg19](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh37_mapping/gencode.v46lift37.basic.annotation.gtf.gz) and [hg38](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz).

### Preprocess Data

1. Run **make_TSS_BED.R** to create a BED file with TSS from the reference genomes, setting the following variables:
   -  **hg19Loc:** The path to the hg19 reference genome
   -  **hg38Loc:** The path to the hg38 reference genome
   -  **hg19BedLoc:** The path to the file where you wish to store the hg19 BED file
   -  **hg38BedLoc:** The path to the file where you wish to store the hg38 BED file
2. Download [bigWigToWig](https://www.encodeproject.org/software/bigwigtowig/).
3. Run the following for each file in BigWig format to convert it to Wig format:
```
./bigWigToWig source.bw dest.wig
```
4. Download [wig2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/wig2bed.html)
5. Run the following for each file in Wig format to convert it to BED format:
```
wig2bed < source.wig | bedmap --echo --max-element -> dest.bed
```
6. Run **fixFormat.R** to fix the format of two BED files, setting the following variables:
   -  **GSM1869138_BJ_PolII_MeDiChi-seq_peaks:** The path to the file containing the sample with this name
   -  **GSM1869138_BJ_PolII_MeDiChi-seq_peaks_out:** The path to the file where you wish to store the version of this sample with the fixed format
   -  **GSM1869150_BJELM_PolII_MeDiChi-seq_peaks:** The path to the file containing the sample with this name
   -  **GSM1869150_BJELM_PolII_MeDiChi-seq_peaks_out:** The path to the file where you wish to store the version of this sample with the fixed format
8. Install bedtools. To do this on a Mac, run the following commands:
```
brew tap homebrew/science
brew install bedtools
```
9. Run the following to find the promoter binding sites for each file:
```
bedtools intersect -wo -a promoters.bed -b peaks.bed > promoter_binding_sites.bed
```

### Generate ChIP-seq Networks

Run **build_ChIPseq_networks.R**, setting the following variables:
-  **baseRepo:** The path to the parent directory where the fibroblast and B lymphocyte BED files are stored


### Evaluate Networks Using FERRET

Run **evaluate_ChIPseq.R**, setting the following variables:
-  **networkDir:** Path to directory where the networks are stored
-  **resultDir:** Path to directory where you wish to store the FERRET results
-  **pathwayFile:** Path to the pathway file downloaded from MSigDB
-  **differentialPathwaysDir:** Path to the directory where you wish to save the differential pathways

## Random Network Evaluation

1. Run **run_FERRET_on_random_networks.R**, setting the following variables:
    -  **inputDir:** The path to the directory where you wish to store the random networks
    -  **outputDir:** The path to the directory where you wish to store the FERRET results
2. Run **consolidate_AUC_random.R**, setting the following variables:
    -  **outputDir:** The path to the directory where the FERRET results are stored

## CPTAC FERRET Evaluation

### Running GRN Inference Methods

#### GRISLI

1. Run **transpose_input.R** to transpose the data, setting the following variables:
    -  **split_dir:** The path to the directory containing the splits
    -  **transpose_split_dir:** The path to where you wish to store the transposed split files
2. Run `./run_grisli_parallel.sh`, which calls **run_grisli_single.m**, setting the following variables:
    -  **EXPRESSION_DIR:** The path to the directory containing the transposed expression data
    -  **PSEUDOTIME_DIR:** The path to the directory containing the pseudotime data
    -  **OUTPUT_DIR:** The path to the directory where you want the output GRNs to be saved
    -  **LOGGING_DIR:** The path to the directory where you want the logs to be saved
    -  **GRISLI_LOCAL:** The path to the directory where GRISLI is installed

#### PIDC
Run `./run_PIDC_parallel.sh`, which calls **runPIDCsingle.jl**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory containing the expression data
-  **OUTPUT_DIR:** The path to the directory where you want the output GRNs to be saved
-  **LOGGING_DIR:** The path to the directory where you want the logs to be saved

#### scSGL

Run `./run_scsgl_parallel.sh`, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory containing the expression data
-  **OUTPUT_DIR:** The path to the directory where you want the output GRNs to be saved
-  **LOGGING_DIR:** The path to the directory where you want the logs to be saved

#### scGeneRai
Run `./run_scgenerai_parallel.sh`, which calls **scgenerai_run_single.py**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory containing the expression data
-  **OUTPUT_DIR:** The path to the directory where you want the output GRNs to be saved
-  **LOGGING_DIR:** The path to the directory where you want the logs to be saved

#### LEAP

Run `./run_LEAP_parallel.sh`, which calls **run_LEAP_single.py**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory containing the transposed expression data
-  **PSEUDOTIME_DIR:** The path to the directory containing the pseudotime data
-  **OUTPUT_DIR:** The path to the directory where you want the output GRNs to be saved
-  **LOGGING_DIR:** The path to the directory where you want the logs to be saved

#### SINGE

1. Run **make_mat_files_singe.R**, setting the following variables:
    -  **expressionPath:** The path to the transposed expression data
    -  **pseudotimePath:** The path to the pseudotime data
    -  **outputPath:** The path where you wish to save the formatted data
2. Run `./run_singe_parallel.sh`, setting the following variables:
    -  **EXPRESSION_DIR:** The path to the directory containing the expression data formatted for SINGE
    -  **SINGE_DIR:** The path to the directory where SINGE is installed
3. Run **convert_mat_to_csv.R**, setting the following variables:
    -  **matPath:** The path where SINGE is installed
    -  **outputPath:** The path where you wish to save the CSV files


### Evaluating Performance Using FERRET
1. Run **format_PIDC_results_for_FERRET.R**, setting the following variables:
    -  **originalDir:** The directory containing the original PIDC results
    -  **consolidatedDir:** The directory containing the reformatted PIDC results
2. Run **run_FERRET.R**, setting the following variables:
    -  **scGeneRai_results:** The path to the directory containing networks inferred using scGeneRai
    -  **scGeneRai_results_formatted:** The path to the directory where you wish to store the reformatted scGeneRai results
    -  **GRISLI_results:** The path to the directory containing the networks inferred using GRISLI
    -  **scSGL_results:** The path to the directory containing the networks inferred using scSGL
    -  **PIDC_results:** The path to the directory containing the networks inferred using PIDC
    -  **LEAP_results:** The path to the directory containing the networks inferred using LEAP
    -  **SINGE_results:** The path to the directory containing the networks inferred using SINGE
    -  **roc_auc_parent_dir:** The path to the parent directory where the FERRET results for each inference method will be stored

### Evaluating Pathways

1. Run **run_FERRET_pathways.R**, setting the following variables:
    -  **scGeneRai_results:** The path to the directory containing networks inferred using scGeneRai
    -  **scGeneRai_results_formatted:** The path to the directory where you wish to store the reformatted scGeneRai results
    -  **GRISLI_results:** The path to the directory containing the networks inferred using GRISLI
    -  **scSGL_results:** The path to the directory containing the networks inferred using scSGL
    -  **PIDC_results:** The path to the directory containing the networks inferred using PIDC
    -  **LEAP_results:** The path to the directory containing the networks inferred using LEAP
    -  **SINGE_results:** The path to the directory containing the networks inferred using SINGE
    -  **pathway_dir:** The path to the parent directory where the pathway analysis results should be stored
    -  **GMT:** The path to the GMT pathway file
2. Run **pathway_frequency_CPTAC_cellTypePairs.R**, to compute the frequency of pathway enrichment when compared across two cell types, setting the following variables:
    -  **pathway_dir_GRISLI:** The path to the directory containing the pathway results for GRISLI
    -  **pathway_dir_LEAP:** The path to the directory containing the pathway results for LEAP
    -  **pathway_dir_PIDC:** The path to the directory containing the pathway results for PIDC
    -  **pathway_dir_scGeneRai:** The path to the directory containing the pathway results for scGeneRai
    -  **pathway_dir_scSGL:** The path to the directory containing the pathway results for scSGL
    -  **pathway_dir_SINGE:** The path to the directory containing the pathway results for SINGE
    -  **pathway_dir_GRISLI_all:** The path to the directory containing the consolidated pathway results for GRISLI
    -  **pathway_dir_LEAP_all:** The path to the directory containing the consolidated pathway results for LEAP
    -  **pathway_dir_PIDC_all:** The path to the directory containing the consolidated pathway results for PIDC
    -  **pathway_dir_scGeneRai_all:** The path to the directory containing the consolidated pathway results for scGeneRai
    -  **pathway_dir_scSGL_all:** The path to the directory containing the consolidated pathway results for scSGL
    -  **pathway_dir_SINGE_all:** The path to the directory containing the consolidated pathway results for SINGE
    -  **pathway_dir_GRISLI_table:** The path to the directory containing the consolidated pathway table for GRISLI
    -  **pathway_dir_LEAP_table:** The path to the directory containing the consolidated pathway table for LEAP
    -  **pathway_dir_PIDC_table:** The path to the directory containing the consolidated pathway table for PIDC
    -  **pathway_dir_scGeneRai_table:** The path to the directory containing the consolidated pathway table for scGeneRai
    -  **pathway_dir_scSGL_table:** The path to the directory containing the consolidated pathway table for scSGL
    -  **pathway_dir_SINGE_table:** The path to the directory containing the consolidated pathway table for SINGE


### Making Bar Plots

Run **consolidate_AUC_by_cell_type_CPTAC.R**, changing the following variables:
-  **rocAucDirScSGL:** Path to the directory storing FERRET results for scSGL
-  **rocAucDirPIDC:** Path to the directory storing FERRET results for PIDC
-  **rocAucDirLEAP:** Path to the directory storing FERRET results for LEAP
-  **rocAucDirScGeneRai:** Path to the directory storing FERRET results for scGeneRai
-  **rocAucDirSINGE:** Path to the directory storing FERRET results for SINGE
-  **rocAucDirGRISLI:** Path to the directory storing FERRET results for GRISLI
-  **ranges_all_file:** Path to the file where you wish to store AUC ranges
-  **jaccard_plot:** Path to the PDF file where you wish to plot Jaccard bars
-  **indegree_plot:** Path to the PDF file where you wish to plot In-Degree bars
-  **outdegree_plot:** Path to the PDF file where you wish to plot Out-Degree bars


## HTAN FERRET Evaluation

### Running GRN Inference Methods

#### SCORPION

1. Download Grand [TF binding](https://granddb.s3.amazonaws.com/tissues/motif/tissues_motif.txt) and [PPI](https://granddb.s3.amazonaws.com/tissues/ppi/tissues_ppi.txt) files.
2. Modify the following **run_scorpion.R** variables:
    -  **motif_file:** Path to the TF binding motif file
    -  **ppi_file:** Path to the PPI file
3. Run **./run_scorpion_HTAN.sh**, setting the following variables:
    -  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
    -  **OUTPUT_DIR:** The path to the directory where you wish to store the output
    -  **LOGGING_DIR:** The path to the directory where you wish to store the logs

#### SCENIC

1. In R, launch the following commands to install SCENIC and save the SCENIC RMD file:
```
if(!require("SCENIC")){  devtools::install_github("aertslab/SCENIC") }
vignetteFile <- file.path(system.file('doc', package='SCENIC'), "SCENIC_Running.Rmd")file.copy(vignetteFile, "SCENIC_myRun.Rmd")
```
2. Download the [motif database](https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/).
3. Modify the following variables in **run_scenic.R**:
    -  **cisTarget_database_dir:** Path to directory where you downloaded the motif database.
4. Run **./run_scenic_HTAN.sh**, setting the following variables:
    -  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
    -  **OUTPUT_DIR:** The path to the directory where you wish to store the output

### Evaluating Performance Using FERRET

Run **run_FERRET_HTAN_AWS.R**, setting the following variables:
-  **networks_SCORPION_file:** The path to the directory containing the SCORPION networks
-  **networks_SCENIC_file:** The path to the directory containing the SCENIC networks

### Evaluating Pathways

1. Run **make_pathways_plot_HTAN.R**, setting the following variables:
    -  **GMT:** The path to the GMT pathway file
    -  **network_dir_SCORPION:** The path to the directory that has the SCORPION networks
    -  **network_dir_SCENIC:** The path to the directory that has the SCENIC networks
    -  **pathway_dir_SCORPION:** The path to the directory where you wish to store the SCORPION pathways
    -  **pathway_dir_SCENIC:** The path to the directory where you wish to store the SCENIC pathways
2. Run **pathway_frequency_HTAN_cellTypePairs.R**, to compute the frequency of pathway enrichment when compared across two cell types, setting the following variables:
    -  **pathway_dir_SCORPION:** The path to the directory containing the pathway results for SCORPION
    -  **pathway_dir_SCENIC:** The path to the directory containing the pathway results for SCENIC
    -  **pathway_dir_SCORPION_all:** The path to the directory containing the consolidated pathway results for SCORPION
    -  **pathway_dir_SCENIC_all:** The path to the directory containing the consolidated pathway results for SCENIC
    -  **pathway_dir_SCORPION_table:** The path to the directory containing the consolidated pathway table for SCORPION
    -  **pathway_dir_SCENIC_table:** The path to the directory containing the consolidated pathway table for SCENIC

### Making Bar Plots

Run **consolidate_AUC_by_cell_type_HTAN.R**, changing the following variables:
-  **rocAucDirSCORPION:** Path to the directory storing FERRET results for SCORPION
-  **rocAucDirSCENIC:** Path to the directory storing FERRET results for SCENIC
-  **ranges_all_file:** Path to the file where you wish to store AUC ranges
-  **jaccard_plot:** Path to the PDF file where you wish to plot Jaccard bars
-  **indegree_plot:** Path to the PDF file where you wish to plot In-Degree bars
-  **outdegree_plot:** Path to the PDF file where you wish to plot Out-Degree bars

## Running Time Evaluation for HTAN

### GRISLI

Run **./run_grisli_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **PSEUDOTIME_DIR:** The path to the directory where the pseudotime data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output
-  **GRISLI_DIR:** The path where GRISLI is stored

### LEAP

Run **./run_LEAP_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **PSEUDOTIME_DIR:** The path to the directory where the pseudotime data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output

### PIDC

Run **./run_PIDC_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output

### SCENIC

Run **./run_SCENIC_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output

### scGeneRai

Run **./run_scgenerai_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output

### SCORPION

Run **./run_scorpion_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output

### scSGL

Run **./run_scsgl_testsingle_HTAN.sh**, setting the following variables:
-  **EXPRESSION_DIR:** The path to the directory where the expression data are stored
-  **OUTPUT_DIR:** The path where you wish to store the output