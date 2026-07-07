# H. pylori Methylome Dynamics Across pH Conditions

Code and data-processing pipeline supporting a study of the DNA methylome of *Helicobacter pylori* 26695 mutant strains grown under neutral and acidic conditions.

## Abstract

*Helicobacter pylori* possess an unusually high number of restriction-modification (R-M) systems relative to its small genome, contributing to a methylome increasingly implicated in bacterial gene regulation. In this study, we analyzed the methylomes of two mutant strains of *H. pylori* 26695: ∆rdxA (control) and ∆rdxA/∆arsS. Each mutant was cultivated under neutral (pH 7) and acidic (pH 5) growth conditions. We identified one conspicuous hypomethylated region of 21 kBp possessing 21 annotated genes across each methylome. Notably, over 600 protein coding regions and 10 different promoters displayed differential methylation between pH conditions, including several virulence factors. The *vacA* gene, encoding the Vacuolating Cytotoxin A, exhibited eight differentially methylated positions between pH 7 and pH 5 within the *H. pylori* 26695 control mutant methylome, potentially contributing to its previously documented 32-fold downregulation of mRNA in acidic environments. pH-dependent methylation changes were widespread within the *cag* pathogenicity island, genes encoding cell envelope proteins including adhesin-encoding *sabA*, *babA*, and *hopQ*, as well as numerous flagellar-associated genes. These results reveal the plasticity of the *H. pylori* methylome and suggest that DNA methylation is responsive to environmental pH in both ArsRS-dependent and independent manners. Methylome dynamics may serve as an important layer of gene regulation in acclimation to hostile gastric environments and promote persistent infection.

## Repository Structure

```
.
├── R files/       # R scripts for statistical analysis and figure generation
├── data/          # Raw and processed sequencing/methylation data
├── deprecated/    # Legacy scripts and files retained for reference, no longer in active use
├── notebooks/     # Jupyter notebooks for exploratory analysis and visualization
├── src/           # Core Python source code (pipelines, utilities, modules)
├── .gitignore     # Files and directories excluded from version control
├── Pipfile        # Python dependency and environment specification (Pipenv)
└── ReadMe.md      # Project documentation (this file)
```

> **Note:** Folder contents above are described at a high level — feel free to edit these descriptions to more precisely match what's inside each directory.

## Getting Started

### Prerequisites

- Python 3.x
- [Pipenv](https://pipenv.pypa.io/) for dependency management
- R (for scripts in `R files/`)

### Installation

Clone the repository and install Python dependencies:

```bash
git clone <[repository-url](https://github.com/skpatterson03/forsyth_lab/tree/main)>
cd <forsyth_lab>
pipenv install
pipenv shell
```

### Usage

1. Place raw methylation/sequencing data in the `data/` directory (see that folder for expected format).
2. Run core processing scripts in `src/` to generate methylation calls and differential methylation results.
3. Use notebooks in `notebooks/` to explore results and reproduce figures.
4. R scripts in `R files/` handle downstream statistical analysis and plotting.

## Data

This project analyzes methylome data from *H. pylori* 26695 ∆rdxA (control) and ∆rdxA/∆arsS mutant strains, each grown under pH 7 (neutral) and pH 5 (acidic) conditions.

## Citation

If you use this code or data, please cite:

> Methylation Dynamics in Helicobacter pylori: Exploring Acidic Stress Effects on Epigenetic Acclimation

## Contact

Sarah Patterson at skpatterson@wm.edu
Dr. Mark Forsyth at mhfors@wm.edu 
