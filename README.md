# Perturb-CITE-seq

Code used in the computational pipeline relating to "Multi-modal pooled Perturb-CITE-Seq screens in patient models define mechanisms of cancer immune evasion"

For any issues, please email cfrangie@broadinstitute.org

If you use this code in your work, please cite: Frangieh, C.J., Melms, J.C., Thakore, P.I. et al. Multimodal pooled Perturb-CITE-seq screens in patient models define mechanisms of cancer immune evasion. Nat Genet 53, 332–341 (2021). https://doi.org/10.1038/s41588-021-00779-1

Count matrices can be downloaded from the Broad Institute Single Cell Portal (SCP1064): https://singlecell.broadinstitute.org/single_cell/study/SCP1064/multi-modal-pooled-perturb-cite-seq-screens-in-patient-models-define-novel-mechanisms-of-cancer-immune-evasion

Raw data can be downloaded from the Broad Institute DUOS system (DUOS-000124): https://duos.broadinstitute.org/

## Overview of Experimental Data

Experimental data was generated according to the following workflow:

![Experimental Data](https://github.com/klarman-cell-observatory/Perturb-CITE-seq/blob/main/experimental_data.png)

## Overview of Computational Workflow

The CRISPR enrichment screen was processed using the [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) software package [1].

The Perturb-CITE-seq data was analyzed according to the following pipeline:

![Computational Pipeline](https://github.com/klarman-cell-observatory/Perturb-CITE-seq/blob/main/computational_pipeline.png)

### Pre-processing of single-cell data

Expression matrices, representing Unique Molecular Identifier (UMI) counts for both scRNA-seq and CITE-seq data were obtained using the Cumulus version 0.14.0 implementation of the CellRanger v.3 workflow with genome reference GRCh38 v3.0.0 and default parameters [2]. Cells with fewer than 200 detected genes or with >18% of detected genes labeled mitochondrial were removed from subsequent analysis. Genes detected in fewer than 200 cells were also removed from further analysis.

## Citations

[1] Li, et al. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology 15:554 (2014)

[2] Li, B., Gould, J., Yang, Y. et al. Cumulus provides cloud-based data analysis for large-scale single-cell and single-nucleus RNA-seq. Nat Methods 17, 793–798 (2020). https://doi.org/10.1038/s41592-020-0905-x
