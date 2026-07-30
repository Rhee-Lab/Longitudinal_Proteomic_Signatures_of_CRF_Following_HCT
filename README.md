# Longitudinal Proteomic Signatures of Cardiorespiratory Fitness Following Hematopoietic Cell Transplantation

 - Associated Publication: (to be filled upon publication)
 - Authors: June-Wha Rhee MD<sup>1</sup>, Lazarina Butkovich PhD<sup>1</sup>, Kavya Achanta PhD<sup>1</sup>, Emma Grigorian<sup>2</sup>, Alysia Bosworth<sup>2</sup>, Sophia Albanese<sup>2</sup>, Alex Flores<sup>2</sup>, Caitlyn Estrada<sup>2</sup>, Mareen Kassabian<sup>2</sup>, Meagan Echevarria MPH<sup>2</sup>, Lanie Lindenfeld MA<sup>2</sup>, Xinyi Du MPH<sup>2</sup>, Congying Xia MD PhD <sup>3,4</sup>, Craig Hyde PhD<sup>5</sup>, Cheryl Keogh-Tow PhD<sup>5</sup>, Cassandra Tierney PhD<sup>5</sup>, Jaron Arbet PhD<sup>2</sup>, Kyuwan Lee PhD<sup>6</sup>, Faizi Jamal MD<sup>1</sup>, F. Lennie Wong PhD<sup>2</sup>, Ryotaro Nakamura MD<sup>8</sup>, James Januzzi, MD<sup>7</sup>, Yun Rose Li<sup>9</sup>, Lee Jones PhD<sup>2</sup>, Bonnie Ky MD, MSCE<sup>3,4</sup>, Raja Mangipudy PhD<sup>5</sup>, Saro Armenian DO MPH<sup>2,x</sup>, Vishal Vaidya PhD <sup>5,x<sup> 
 - Affiliations:
    - <sup>1</sup> Department of Medicine, City of Hope Comprehensive Cancer Center, Duarte, CA 
    - <sup>2</sup> Department of Population Sciences, Beckman Research Institute of the City of Hope, Duarte, CA
    - <sup>3</sup> Division of Cardiology, Department of Medicine, Perelman School of Medicine, University of Pennsylvania, Philadelphia, PA 
    - <sup>4</sup> Thalheimer Center for Cardio-Oncology, Abramson Cancer Center, Perelman School of Medicine at the University of Pennsylvania, Philadelphia, PA 
    - <sup>5</sup> Pfizer Research and Development, 1 Portland Street, Cambridge, MA 
    - <sup>6</sup> Department of Kinesiology and Sports Studies, Ewha Womans University, Seoul, South Korea 
    - <sup>7</sup> Division of Cardiology, Department of Medicine, Massachusetts General Hospital, Boston 
    - <sup>8</sup> Department of Hematology & Hematopoietic Cell Transplantation, City of Hope Comprehensive Cancer Center, Duarte, CA 
    - <sup>9</sup> Department of Radiation Oncology, City of Hope Comprehensive Cancer Center, Duarte, CA 
    - <sup>x</sup> Co-Corresponding Authors

 - GitHub repository author: Lazarina Butkovich, PhD

## Study Overview
(to be filled upon publication)

## Repository Setup
- This respository uses [renv](https://rstudio.github.io/renv/articles/renv.html) as the R package manager. Package versions are recorded in the renv.lock file, so the project environment can be reproduced on another machine (Windows).
- Prerequisites:
    - Install a recent version of R (used R 4.5.2).
- Clone the repository:
```
git clone https://github.com/lbutkovich/Longitudinal_Proteomic_Signatures_of_CRF_Following_HCT.git
cd <Longitudinal_Proteomic_Signatures_of_CRF_Following_HCT>
```
- Install renv:
```
install.packages("renv")
```
- Restore the project environment using the renv.lock file
```
renv::restore()
```
- Verify installation:
```
renv::status()
```

## Running the Workflow
- To run all proteomics workflow: run_all_proteomics_analysis.R
- To run individual scripts, see the additional .R files in the proteomics_scripts folder.
- To change workflow parameters, see config.R
- Script details:
    - (1) proteomics_batch_correction_ComBat.R
        - combine Olilnk NPX data across multiple plates
        - perform batch correction across plates, using the ComBat tool
        - generate mean NPX for PTID-visit replicates
    - (2) proteomics_data_processing.R
        - Plot normalized protein expression (NPX) distribution, pre- and post-filters
        - Plot variance distribution
        - Apply filters:
            - High missingness
            - Low-variance
        - Calculate delta VO2 values (6m-Baseline and 12m-Baseline VO2)
        - Export processed data and patient metadata
        - Plot PCA
        - Plot variance distribution by timepoint
    
    Proteomics Linear Regression Approach:
    - (3) proteomics_linear_regression.R
        - Fit linear regression model per contrast
        - Generate p-value histograms
        - Plot beta coefficient distributions
    - (4) proteomics_linear_regression_enrichment.R
        - Perform GSEA and ORA per contrast
    - (5) proteomics_linear_regression_plots.R
        - Generate GSEA bar plots
        - Generate time series heatmaps (protein-specific)
        - Generate time series heatmaps (for annotations)
        - Generate volcano plot of proteins

    Proteomics Top vs. Bottom Quartile VO2, Linear Regression Approach:
    - (6) proteomics_top_vs_bottom_diff_exp.R
        - Rank patients by VO2peak per timepoint; split into top and bottom 25%
        - Differential expression analysis between extreme quartiles using limma
        - Generate p-value histograms and logFC distributions
        - Export results as Excel files
    - (7) proteomics_top_vs_bottom_diff_exp_enrichment.R
        - Perform GSEA and ORA per contrast
    - (8) proteomics_top_vs_bottom_diff_exp_plots.R
        - Generate GSEA bar plots
        - Generate volcano plots
        - Generate time series heatmaps (protein-specific and annotation-level)

    Additional files:
    - generate_additional_figures.R
        - create additional and ad hoc plots
    - common_functions.R
        - Compilation of functions shared across multiple scripts
    - config.R
        - Location of values to set manually

- Input files:
    - Contact Lazarina Butkovich (lbutkovich@coh.org) to inquire about access to raw data for input to this workflow.