# **Precision probabilistic mapping reveals personalized functional topography in the early postnatal human brain**

E-mail: [zhaojianlong@mail.bnu.edu.cn](mailto:zhaojianlong@mail.bnu.edu.cn); [tengdazhao@bnu.edu.cn](mailto:tengdazhao@bnu.edu.cn)

This repository provides code and source data that support the findings of the article entitled "Precision probabilistic mapping reveals personalized functional topography in the early postnatal human brain“ 

## Overview

This repository includes standalone software, source code, and demonstration data used in the analyses.
 All data required to reproduce the findings are publicly available, including:

- Group- and individual-level functional topographies
- Individual variation in neonatal brain networks
- Functional topography undergoes refinement in early postnatal period and encodes brain maturity
- Functional topographic refinement in neonatal brain is anatomically constrained by the structural maturation 
- Altered network topography and accelerated functional maturation in preterm infants
- Source data for figure visualization

All files are hosted in a publicly accessible cloud repository (https://github.com/zhaohuaxishi1/neonatal_Individual_Parcellation).

### **Public datasets used**

- **dHCP dataset:** https://nda.nih.gov/
- **HCP dataset:** https://www.humanconnectome.org/

Source data supporting the findings are provided with this paper.

### Data

All data are available in Releases (https://github.com/zhaohuaxishi1/neonatal_Individual_Parcellation/releases/tag/v1.0.0).


- `group_icn_soft`, `group_icn_hard/`: Group-level soft and hard functional topographies
- `individual_parcellation/`: Individualized functional parcellations 
- `individual_variability/`: Voxelwise and network-level individual variability
- `network_loading/`: Feature of Age & behavioral outcome scores at 18 months
- `structure_basis/`: Structural features & metrics (`sub_metric_all_term.mat`)
- `classification/`: Matched term-preterm data and labels (`V_vector_FDAll.mat`, `ClassificationlabelFDAll.mat`)
- `*.mat` (e.g., `PredictionAgeScore.mat`): Age & behavioral outcome scores at 18 months
- `init.mat`:  soft functional topographies (Voxel, Network number)
- `Resliced3mm_GmMask_final.nii`: Final gray matter mask

### Code (`/code`)

| **Module** | **Functionality**                                          | **Scripts**                                                  |
| ---------- | ---------------------------------------------------------- | ------------------------------------------------------------ |
| S1         | Generate group and individual functional topographies      | `s1_1_Group_level_probabilistic_network_topographies.m`, `s1_2_Hard_parcellation.m` |
| S2         | Analyze individual variability and boundary concentraction | `s2_1_Voxel_wise_individual_variability.m`, `s2_2_Network_level_individual_variability.m`, `s2_3_1_Boundary_concentration.m`, `s2_3_2_Plot_Boundary_concentration.R` |
| S3         | Predict brain age using topographic features               | `S3_1_AgeEffect_totalLoading.R`, `S3_2_AgeEffect_voxelwise_networks.R`, `S3_3_Predict_brain_age.m` |
| S4         | Predict cognitive, language, and motor outcomes            | `S4_1_Predict_cognitive_score.m`, `S4_2_Predict_language_score.m`, `S4_3_Predict_motor_score.m` , `S4_4_Predict_qchat_score.m` |
| S5         | Structural basis of functional development                 | `S5_1_GAM_structural_maturation.R`, `S5_2_Visualize_structural_age_effects.m`, `S5_3_Struc_Func_Coupling.m` , `S5_4_network_level_PLS_All.m`|
| S6         | Classify preterm vs term infants & visualize ROC           | `S6_1_Classify_preterm_vs_term.m`, `S6_2_Visualize_classification_ROC.m`, `S6_3_Visualize_preterm_BAG_violin.m` |

### Key Analyses & Workflow

#### 1. **Group and Individual Functional Topographies**

- **Group-level networks**: computed via NMF on the dHCP dataset
- **Individual parcellation**: subject-specific parcellations using a NMF

#### 2. **Individual Variability**

- **Voxelwise and network-level variability** quantified across subjects
- **Boundary concentration** maps identify topographic transition zones

#### 3. **Brain Maturity Prediction**

- **Age prediction** via Ti-PCA + SVR on network-level features
- **Voxelwise GAM** used to assess spatial maturation patterns

#### 4. **Neurodevelopmental Outcome Prediction**

- **18-month outcomes** (cognitive, language, motor, Q-chat) predicted using functional topographical features

#### 5. **Structural Basis of Functional Development**

- Cortical metrics (myelination, sulcal depth, curvature, thickness) modeled with GAM
- Significant maturational regions visualized via GIFTI surfaces
- partial least-square (PLS) analyses to assess the latent relationship between these maturation maps (F values of the age effect) and the contribution weights of network-specific functional topography.
#### 6. **Preterm Functional Alterations**
- **Classification (term vs preterm)** using SVM with Ti-PCA features
- ROC curves and decision values are visualized to assess model performance
- Personalized brain age gap (BAG) analysis
# Software & Tools

All analyses were conducted using open-source software and publicly available toolboxes.

### **Preprocessing**

- dHCP functional pipeline: https://git.fmrib.ox.ac.uk/seanf/dhcp-neonatal-fmri-pipeline/-/tree/master
- dHCP structural pipeline: https://github.com/BioMedIA/dhcp-structural-pipeline
- HCP minimal preprocessing pipeline: https://github.com/Washington-University/HCPpipelines/releases/

#### Postprocessing

- SPM12 toolbox: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- GRETNA toolbox v2.0.0: https://www.nitrc.org/projects/gretna/
- CIFTI-Matlab toolbox v2: https://github.com/Washington-University/cifti-matlab/
- MATLAB R2020b: https://www.mathworks.com/products/matlab.html

#### Analysis

- Individual Parcellation method: https://github.com/hmlicas/Collaborative_Brain_Decomposition
- Support Vector Regression (SVR) model: https://github.com/ZaixuCui/Pattern_Regression_Clean and LIBSVM (v3.25): https://www.csie.ntu.edu.tw/~cjlin/libsvm/
- Support Vector Classification: https://github.com/ZaixuCui/Pattern_Classification
- Generalized Additive Models (GAM): implemented in R using the mgcv package (https://cran.r-project.org/web/packages/mgcv/index.html)
- Partial Least Squares (PLS): https://github.com/valkebets/myPLS-1

#### **Visualization**

- BrainNet Viewer: [www.nitrc.org/projects/bnv](http://www.nitrc.org/projects/bnv))
- Connectome Workbench: https://www.humanconnectome.org/software/connectome-workbench-R 
- ggplot2: https://ggplot2.tidyverse.org/
- GIFTI and MATLAB (R2020b)

## System Requirements

#### Linux Server Cluster

- Ubuntu 22.04 LTS
- Intel® Xeon® CPU E5-2683 v4 @ 2.10 GHz
- 128 GB RAM
- 32 CPU cores
- 24 compute nodes

####  Windows Workstation

- Microsoft Windows 10 Enterprise (64-bit)
- Intel® Core™ i7-9700 CPU
- 24 GB RAM
- 8 CPU cores

#### Runtime

The complete analysis pipeline was run on a 24-node Linux computing cluster and required approximately 15 days to complete all major experiments and generate the primary results reported in the manuscript.

## 🧠 Output Examples

| **Analysis**                        | **Output Folder**                                            |
| ----------------------------------- | ------------------------------------------------------------ |
| Group-level soft atlas              | `output/fig_grp_soft/`                                       |
| Individual variability & boundaries | `output/individual_variability/`, `output/variability_loading_corr/` |
| Brain age prediction                | `output/gam_totalloading/`, `output/gam_voxelwise/`          |
| Cognitive prediction                | `output/SVR_10CV_TiPCA_Cognitive/`                           |
| Structural maturation               | `output/gam_structurematuration/`                            |
| Classification ROC                  | `output/classification/SVM_10CV_TiPCA/`                      |

