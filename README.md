# LABM Multi-Camera Peak Analysis Pipeline

Automated image analysis pipeline for the **Large Angle Beamstrahlung Monitor (LABM)** 8-camera system — part of the [Belle II experiment](https://www.belle2.org/) at the SuperKEKB accelerator in Tsukuba, Japan.

This pipeline processes images from all 8 LABM cameras across a full accelerator run, automatically detects photon peaks in every image, identifies the highest-intensity events per camera, generates diagnostic heatmaps, and exports a structured ROOT ntuple linking camera data to accelerator beam conditions for downstream analysis.

---

## What it does

Each camera captures an image of the photon distribution on its detector face during accelerator operation. For a full run, this can mean tens of thousands of images across 8 cameras, each associated with a specific set of beam conditions. This pipeline automates the entire analysis chain:

1. **Reads accelerator condition data** — ingests timestamped beam parameters for a full run: HER/LER currents, luminosity (HER, LER, ECL), and beam position monitor readings across 8 magnet locations in both X and Y for each ring (22 parameters per event)

2. **Processes every image, every camera** — for each of the 8 cameras, converts the corresponding PNG image into a 2D intensity heatmap from pixel ARGB values, applies 60×60 rebinning for noise reduction, and runs ROOT's `TSpectrum2` peak-finding algorithm to detect up to 400 photon peaks per image

3. **Tracks peak statistics per camera, per event** — stores the number of peaks found and the position/intensity of each peak for every camera across every event in the run

4. **Identifies the top events automatically** — for each of the 8 cameras, finds the 5 events with the highest peak intensity across the entire run, without manual inspection

5. **Generates targeted diagnostic heatmaps** — produces full heatmap PDFs (contour + 3D surface, with detected peaks marked) only for those top-intensity events, keeping output focused on the most informative frames instead of generating thousands of redundant plots

6. **Builds a full diagnostic histogram set** — produces 1D histograms for all 22 accelerator parameters across the run (currents, luminosities, magnet positions) for beam condition analysis

7. **Exports a ROOT ntuple** — writes a structured tree (`labmtree`) with one entry per event, containing all beam condition parameters, the image filename for each of the 8 cameras, and the peak count detected per camera — ready for further analysis in ROOT, Python (uproot), or any ROOT-compatible tool

---

## Example output

- **Heatmap PDFs** (top 5 events per camera): two-page PDF per event showing a 2D contour map with peak markers and a 3D surface plot
- **Histogram PDFs**: full set of beam condition distributions for the run
- **ROOT ntuple** (`run6.root`): structured tree with 30 branches — 21 beam condition floats, 1 timestamp, 8 camera image filenames, and 8 per-camera peak counts — one row per event

---

## File structure

```
LABM-Camera-Analysis/
├── peaks.C                        ← Main script: full pipeline (V3)
├── README.md
└── runs/                          ← Run directory (maindir)
    └── run6/                      ← Run 6 directory (example)
        ├── data-run6-all-...      ← Accelerator condition data file
        ├── pngPics/               ← Input camera images (8 cameras)
        ├── images/                ← Output: heatmap PDFs for top events
        ├── graphs/                ← Output: beam condition histogram PDFs
        └── ntuple/
            └── run6.root          ← Output: structured ROOT ntuple
```

---

## Requirements

- [ROOT](https://root.cern/) 6.x or later
- ROOT compiled with image processing support (`TASImage`)
- Camera PNG images from all 8 LABM cameras
- Accelerator condition data file (structured text format, 22 parameters per event)

---

## How to run

```bash
root -l -b -q peaks.C
```

- `-l` skips the ROOT splash screen
- `-b` runs in batch mode (no graphical windows)
- `-q` closes ROOT when finished

Before running, update the directory paths near the top of `peaks.C`:

```cpp
// Main run directory
string maindir = "/path/to/your/run/";

// Camera PNG images directory (relative to maindir)
string pngdir = maindir + string("pngPics/");

// Heatmap PDF output directory
string heatmdir = maindir + string("images/");

// Histogram PDF output directory
string histdir = maindir + string("graphs/");

// Accelerator data file
sprintf(n, "/path/to/your/data-file");

// Output ntuple location
string rootfilenamedummy = maindir + string("ntuple/run6.root");
string treenamedummy = string("run6");
```

To control how frequently images are processed during the first pass, adjust the `div` parameter in the `readevent()` call inside `peaks()`. The script supports a maximum of `maxeventnum` (30,000) events and up to `maxpeaks` (400) peaks per image — both adjustable via the `#define` constants at the top of the file.

---

## Context

This pipeline was developed as part of the LABM camera upgrade at KEK, where the original PMT-based readout was expanded to an 8-camera system covering all optical channels of the detector. Automatically scanning every image for every event in a run — and surfacing only the highest-intensity frames per camera — makes it practical to monitor beam-induced photon distributions across a full run without manual review of tens of thousands of images.

The ROOT ntuple output bridges the camera analysis with standard accelerator physics analysis workflows, allowing peak counts and intensities to be directly correlated with beam currents, luminosity, and magnet positions.

---

## Skills demonstrated

`C++` · `ROOT/CERN` · `Image processing` · `2D peak detection` · `Large-scale automated data pipelines` · `Multi-source data integration` · `Structured data export (ROOT ntuples)` · `Scientific visualization` · `Beam diagnostics`

---

## Author

**Daniel Enrique Ricalde Herrmann**  
PhD in Physics — Wayne State University (2024)  
daniel.ricalde.h@gmail.com
