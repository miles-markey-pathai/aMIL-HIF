## pathai-aMIL+HIF: Data-driven interpretability analysis of multiple instance learning models

# Overview
pathai-aMIL+HIF is the primary code repository for reproducing analyses in: "Spatial mapping of immunosuppressive CAF gene signatures in H&E-stained images using additive multiple instance learning"

# Installation
1. Clone this repo to your local machine using https://github.com/Path-AI/aMIL-HIF.git
2. Estimated install time: <1 minute
3. Estimated run time: varies between scripts (1-5 minutes)
# Contents
* The main directory contains all the notebooks and utility python files to reproduce the analysis in this manuscript
* /data contains all raw and cached data objects needed to reproduce the analysis, as well as to read in the training and evaluation data included in the manuscript.
# Version and package requirements
* To install the Python packages and dependencies needed to operate this code, please use Anaconda or Miniconda. From within this directory, do:
* conda env create --name amil-hif --file=environment.yml
* conda activate amil-hif

# Abstract
The relative abundance of cancer-associated fibroblast (CAF) subtypes influences a tumor’s response to treatment, especially immunotherapy. However, the extent to which the underlying tumor composition associates with CAF subtype-specific gene expression is unclear. Here, we describe an interpretable machine learning (ML) approach, additive multiple instance learning (aMIL), to predict bulk gene expression signatures from H&E-stained whole slide images (WSI), focusing on an immunosuppressive LRRC15+ CAF-enriched TGFβ-CAF signature. aMIL models accurately predicted TGFβ-CAF across various cancer types. Tissue regions contributing most highly to slide-level predictions of TGFβ-CAF were evaluated by ML models characterizing spatial distributions of diverse cell and tissue types, stromal subtypes, and nuclear morphology. In breast cancer, regions contributing most to TGFβ-CAF-high predictions (“excitatory”) were localized to cancer stroma with high fibroblast density and mature collagen fibers. Regions contributing most to TGFβ-CAF-low predictions (“inhibitory”) were localized to cancer epithelium and densely inflamed stroma. Fibroblast and lymphocyte nuclear morphology also differed between excitatory and inhibitory regions. Thus, aMIL enables a data-driven link between tissue phenotype and transcription, offering biological interpretability beyond typical black-box models.

# License
The contents of this repository are made available under CC BY-NC 4.0 (as provided in license.txt).


