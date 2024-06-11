# NISARA
## Network Analysis of Structure Similarity Alignment V1.0

NISARA is a simple Python-based script for clustering structures. Instead of all-vs-all BLAST (E-value) comparisons, NISARA uses all-vs-all TM-align comparisons (TM-score). NISARA produces a CLANS file which can be visualized in the CLANS program.

## Installation

**Requirements:** Anaconda installed on the target computer. An version of Anaconda can be downloaded from: https://www.anaconda.com/.

### From Anaconda repository:

1.	Create a clean conda environment: `conda create -n clans_2_0`

2.	Activate the newly created environment: `conda activate clans_2_0`

3.	Install the clans package from Anaconda repository by using the following command:

      **On MacOS / Windows:**

      `conda install -c inbalpaz clans -c defaults -c conda-forge`

      **On Linux:**

      `conda install -c inbalpaz clans_linux -c defaults -c conda-forge`

4. **On Linux only** run the following command: `pip install PyQt5`

### From source using conda:

1. Download CLANS latest release from: https://github.com/inbalpaz/CLANS/releases.

2. Extract the tar.gz file into the desired working-directory.

3. Create a new conda environment using the ‘clans_2_0.yml’ file (located in the root directory of CLANS) using the following command:

   **On Linux:** `conda env create -f NISARA_v1.yaml`

4. Activate the newly created environment: `conda activate clans_2_0`

#### Set the display server on Linux
In case the default display server on your Linux distribution is Wayland, you should switch to Xorg (x11) in order to enable CLANS to work properly.

1. Check the currently used display server: `echo $XDG_SESSION_TYPE`

2. If the result is ‘wayland’, you should switch to ‘x11’. In most popular linux distributions this can be done by uncommenting WaylandEnable=false in the /etc/gdm3/custom.conf file.
    
## Usage

### Open the GUI-based visualisation tool

Within the activated clans_2_0 conda environment, type:

`python -m clans [-load <network file path>] [options]
`

When clans is executed without an input-file, the GUI is opened empty and an input-file can be loaded from the ‘File’ menu.

### Running CLANS in command-line mode

The command-line mode is executed using the ‘-nogui’ flag option. It can be used to perform a BLAST search in order to create a matrix of sequence similarities and/or to perform a specific number of iterations of the force-directed graph layout calculation (which can be later loaded and displayed in the visualizing tool). It is also recommended in cases of large datasets, where the clustering can be done in the background and the resulted clans map can later be loaded and explored using the graphical interface.

Within the activated clans_2_0 conda environment, type:

`python -m clans -nogui -infile <fasta_file_path> -saveto <destination_file_path> [options]
`

or

`python -m clans -nogui -load <network_file_path> -dorounds <number of iterations> -saveto <destination_file_path> [options]`

## Tutorial

A detailed tutorial is found here: https://github.com/inbalpaz/CLANS/tree/master/clans/manual/Manual.pdf

or can be opened from CLANS visualisation tool (Help -> CLANS manual).
