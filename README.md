# Bioinformatics Data Integration: Protein–Metabolite Interaction Analysis

## Overview

During my Bioinformatics Research Internship at Huntsman Cancer Institute, I developed an R-based pipeline to integrate MIDAS experimental data with the Human Metabolome Database (HMDB) and UniProt annotations. The goal was to identify protein–metabolite interactions (PMIs), distinguish previously annotated interactions from those not found in existing database annotations, and explore patterns in the results.

By combining experimental data with biological databases, this project helped identify potential new protein–metabolite associations and characterize interaction patterns for further investigation.

## Key Results

* **86 significant protein–metabolite interactions** identified in the final analysis.
* **80 interactions were previously unannotated** in the databases examined, highlighting candidates for further investigation.
* Identified notable outliers in PCA of interaction profiles, including Propionyl-CoA and HSPA5.
* Developed visualizations to examine interaction significance, frequently interacting proteins and metabolites, and variation in interaction profiles.

## Methods

### 1. Data Integration and Annotation

* Combined positive- and negative-mode MIDAS metabolomics data.
* Matched metabolite identifiers against HMDB records.
* Compared UniProt protein identifiers with HMDB protein associations to classify interactions as known or not previously annotated in the database.
* Used R to clean, filter, and organize interaction-level data.

### 2. Statistical Filtering

Applied significance thresholds to identify interactions for further analysis:

* Log2 fold change > 1
* q-value < 0.01

### 3. Exploratory Data Analysis and Visualization

* **Interaction classification:** Compared previously annotated and unannotated interactions.
* **Frequency analysis:** Identified proteins and metabolites with the highest numbers of significant interactions.
* **Principal component analysis (PCA):** Constructed protein-by-metabolite and metabolite-by-protein matrices using composite interaction scores to explore variation in interaction profiles.
* **Volcano plots:** Visualized fold changes and statistical significance for metabolite interactions with a selected protein.

The composite interaction score was calculated as:

`log2 fold change × −log10(p-value)`

## Technologies

* **Programming:** R
* **Data manipulation:** dplyr, tidyverse, stringr
* **Visualization:** ggplot2, ggrepel
* **Dimensionality reduction:** Principal component analysis (PCA)
* **Biological databases:** UniProt, HMDB
* **Data processing:** Statistical filtering, identifier matching, data transformation, and missing-value handling

## Significance and Limitations

Integrating MIDAS experimental results with existing biological annotations provides a way to investigate protein–metabolite interactions at scale. The large proportion of interactions not found in the databases examined highlights potential gaps in existing annotations.

However, an interaction classified as previously unannotated is not necessarily a confirmed novel biological interaction. Experimental validation and further database comparisons would be needed to establish its biological significance.

PCA provided an exploratory view of interaction-profile variation, while the identified outliers offer potential directions for further investigation.

## Future Directions

Potential extensions include incorporating additional biological resources, such as BRENDA and STITCH, to expand annotation coverage and investigate a broader range of protein–metabolite associations.

## Acknowledgments

This work was completed during my Bioinformatics Research Internship at Huntsman Cancer Institute, University of Utah, with guidance and support from the research team.

**Project author:** Kainoa Martinez
