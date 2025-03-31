
# ProteoPlotter

## Overview
ProteoPlotter is an interactive web-based software developed to complement [Perseus](https://maxquant.net/perseus/) in generating multidimensional figures for the analysis of proteomics datasets.  

## Features
- **Multi-Dimensional Visualizations**: Create [1D annotation](https://doi.org/10.1186/1471-2105-13-S16-S12) heat maps, Venn diagrams, volcano plots, principal component analysis (PCA) plots with data ellipses, UpSet plots, and dynamic range plots.
- **User-Friendly Interface**: ProteoPlotter offers an intuitive graphical user interface built within the R-Shiny framework.
- **PNG Download**: Export your visualizations in a high-quality format.
- **Comprehensive Guide**: Instructions for preparing and uploading Perseus output files are provided to facilitate ease of use.

## Getting Started
1. **Download the `ProteoPlotter_Windows.zip` folder**: The latest version is available [here](https://github.com/JGM-Lab-UoG/Extending_Perseus-Esther-/releases/tag/v1.0.0).
2. **Access the folder**: Unzip the folder, then click the `run.vbs` file to run ProteoPlotter on your PC. 
3. **Prepare data or apply sample datasets**: Generate tabular output files from Perseus or utilize the sample input files provided.
4. **Upload and sisualize**: Upload input files to the ProteoPlotter software to create data visualizations!
5. **Visit the `Guide` tab on ProteoPlotter**: Feel free to peruse instructions under this tab for more help with uploading data and utilizing the software.
6. **Source code**: The source code for ProteoPlotter is available within the `ProteoPlotter_RScript` folder. You can download the folder to run the software through RStudio.

## OS Compatibility
ProteoPlotter is compatible with the Windows OS. 

## Sample Input Datasets
Within the `Sample_Files` folder of this repository, we have provided sample datasets exported from Perseus (v2.1.3.0) for visualizing the dynamic proteome of *Klebsiella pneumoniae*, based on a [previous study](https://doi.org/10.3389/fmicb.2020.00546). The Perseus session file for the analytical workflow that generates the sample datasets is also available, allowing users to follow the initial data analysis process. 
- `Sample_1d_matrix.txt`: Input file for generating 1D annotation enrichment heatmaps.
- `Sample_protein_main.txt`: Input file for generating volcano, PCA, and dynamic range plots.  
- `Sample_curve_matrix.txt`: Input file for overlaying false discovery rate (FDR) threshold curves on volcano plots. 
- `Sample_venn_matrix.txt`: Input file for generating Venn diagrams and UpSet plots. 
- `Sample_file_session_Perseus2.1.3.0.sps`: Perseus data analysis session file which can be opened within the Perseus platform. 



