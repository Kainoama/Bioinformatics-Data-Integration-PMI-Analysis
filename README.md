# Protein–Metabolite Interaction Analysis Through Biological Database Integration

## Overview

During my Bioinformatics Research Internship at Huntsman Cancer Institute, I developed the R-based backend pipeline for a project integrating MIDAS experimental data with annotations from the Human Metabolome Database (HMDB) and UniProt.

MIDAS (Mass Spectrometry Integrated with Equilibrium Dialysis for the Discovery of Allostery Systematically) identifies protein–metabolite binding events. By comparing MIDAS results with existing biological annotations, the pipeline helped distinguish interactions already represented in the databases from those that were not previously annotated, supporting hypothesis generation and follow-up research.

## Key Findings

* Identified **86 significant protein–metabolite interactions** in the final analysis.
* **80 of the 86 interactions were not previously annotated** in the databases examined.
* Used principal component analysis (PCA) to explore variation in protein–metabolite interaction profiles. Propionyl-CoA and HSPA5 appeared as notable outliers in the analyses.
* Examined significant interactions and their distributions through bar charts, PCA plots, and a volcano plot focused on albumin (ALB).

Interactions that were not found in the databases were treated as potentially novel candidates, not as experimentally confirmed new biological interactions.

## Workflow

### 1. Data Preparation and Integration

* Loaded MIDAS statistics and negative-mode metabolomics metadata.
* Matched negative-mode metabolite identifiers against HMDB records.
* Retrieved HMDB protein-association annotations for matched metabolites.
* Compared protein identifiers in the MIDAS results with associated UniProt identifiers to classify interactions as previously annotated or not found in the examined annotations.

### 2. Statistical Filtering

Filtered interactions for downstream analysis using:

* **Log2 fold change > 1**
* **q-value < 0.01**

These thresholds were used to identify significant interactions for the main known-versus-unannotated comparison and related visualizations.

### 3. Exploratory Analysis and Visualization

* **Interaction counts:** Compared previously annotated and unannotated interactions.
* **Frequency analysis:** Identified proteins and metabolites with the highest numbers of significant interactions.
* **PCA:** Built protein-by-metabolite and metabolite-by-protein matrices from composite interaction scores to explore differences in interaction profiles.
* **Volcano plot:** Visualized metabolite interactions with albumin using log2 fold change and statistical significance.

The composite score used in the PCA was calculated as:

`log2 fold change × −log10(p-value)`

PCA was used as an exploratory method to summarize variation in interaction profiles. Clustering and outliers were treated as candidates for further investigation rather than proof of biological function.

## My Contributions

* Developed the R-based backend pipeline for integrating MIDAS data with HMDB and UniProt annotations.
* Implemented metabolite matching, interaction classification, statistical filtering, and downstream analysis.
* Created visualizations to summarize interaction status, interaction frequency, PCA results, and statistical significance.
* Collaborated on a project that also included an interactive Shiny application. The front-end interface was developed by my co-intern; my contribution focused on the backend analysis pipeline.

## Tools and Technologies

* **Language:** R
* **Data wrangling:** `dplyr`, `tidyverse`, `stringr`
* **Visualization:** `ggplot2`, `ggrepel`
* **Dimensionality reduction:** Principal component analysis (PCA)
* **Biological databases:** UniProt, HMDB
* **Methods:** Data integration, identifier matching, statistical filtering, exploratory data analysis

## Future Directions

Potential extensions include incorporating additional resources such as BRENDA and STITCH to broaden annotation coverage and support analysis of larger, more diverse datasets.

## Acknowledgments

This project was completed during my Bioinformatics Research Internship at Huntsman Cancer Institute, University of Utah, with mentorship from Paul Stewart and Kevin Hicks and collaboration with McKenna Wilson.

**Author:** Kainoa Martinez

