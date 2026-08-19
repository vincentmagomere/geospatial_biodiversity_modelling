# Insect-Based Biodiversity Intactness (IBI)
### A reproducible workflow for modelling insect diversity and computing biodiversity intactness across Africa

This repository contains all scripts, workflows, and resources used for the study:

**“An insect- and EO-based indicator reveals microhabitat-driven biodiversity intactness across Africa.”**

The repository provides a complete implementation of the insect-based Biodiversity Intactness (IBI) indicator, including:

- Preprocessing of insect occurrence data  
- Derivation of environmental predictor layers  
- Random Forest species distribution modelling  
- Per-pixel biodiversity intactness computation  
- Uncertainty quantification  
- Visualizations and statistical comparisons  

---

## 📁 Repository Structure


---

## 🔧 Requirements

### Software
- **R ≥ 4.2**
- **Python ≥ 3.9**
- **QGIS ≥ 3.22**
- **Google Earth Engine account**

### Key R packages
- `terra`, `raster`, `sf`, `dismo`
- `spThin`
- `ggplot2`, `dplyr`

### Key Python packages
- `scikit-learn`
- `rasterio`, `geopandas`
- `numpy`, `pandas`

---

## 📘 Project Overview

This repository implements the full workflow for creating the **Insect-Based Biodiversity Intactness (IBI)** indicator.

### 1. Current insect diversity (D)
Modelled using a **Random Forest** classifier trained on:
- Coleoptera, Lepidoptera, Odonata occurrences  
- Sentinel-2 spectral metrics  
- GEDI canopy height  
- SRTM elevation  
- Pseudo-absence/background points  

**Spatial block cross-validation** ensures no spatial leakage.

---

### 2. Potential insect diversity (P)
Baseline values for **pre-human potential diversity**, assigned per habitat class.

Two scenarios included:
- **High P**
- **Low P**

Used for sensitivity analysis.

---

### 3. Human footprint (hF)
Derived from the **Global Human Modification (GHM)** dataset, normalized 0–1.

---

### 4. Biodiversity Intactness Indicator (IBI)

\[
\text{IBI} = \frac{D}{P} \times (1 - hF)
\]

Outputs range from **0 (degraded)** to **1 (intact)**.

---

### 5. Uncertainty
Combined from:
- Distance-to-occurrence  
- Predicted suitability  
- Spatial block CV performance  

Used to identify data-poor or model-uncertain regions.

---

## 🚀 Running the Workflow

### **Step 1 — Prepare Occurrence Data**
```bash
Rscript 1_data_preprocessing/occurrence_cleaning.R
Rscript 1_data_preprocessing/spatial_thinning.R
Rscript 1_data_preprocessing/pseudo_absence_generation.R
