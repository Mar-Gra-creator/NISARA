# NISARA
## Network Analysis of Structure Similarity Alignment V1.0d

NISARA is a simple Python-based script for clustering structures. Instead of all-vs-all BLAST (E-value) comparisons, NISARA uses all-vs-all TM-align comparisons (TM-score). NISARA produces a CLANS file which can be visualized in the CLANS program.

## Installation

**Requirements:** Anaconda installed on the target computer. An version of Anaconda can be downloaded from: https://www.anaconda.com/.

### From source using conda:


1. Copy repo (extract input) into the desired working-directory.

2. Create a new conda environment using the `NISARA_v1.0.yml` file (located in the root directory of NISARA) using the following command:

   **On Linux:** `conda env create -f NISARA_v1.0.yaml`.

3. Activate the newly created environment: `conda activate NISARA_v1.0`.

4. Run `python NISARA_v1.0.py` to display description.

5. Simple run: `python NISARA_v1.0.py nucleases 0.4 0.7 2 8`

## Tutorial

A detailed tutorial is found here: ...


