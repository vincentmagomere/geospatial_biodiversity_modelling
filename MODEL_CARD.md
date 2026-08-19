# Random Forest Insect Diversity Model — Model Card
_Africa-wide ecological suitability (D) model for the Insect-Based Biodiversity Intactness (IBI) Indicator_

**Authors:** Tobias Landmann, Vincent Magomere, Komi Mensah Agboka  
**Affiliation:** International Centre of Insect Physiology and Ecology (icipe)  
**Contact:** tlandmann@icipe.org  
**Repository:** https://github.com/icipe-official/Biodiversity-Project  

---

## 🧭 Model Overview

The Random Forest Insect Diversity Model estimates **current insect ecological suitability (0–1)** across Africa and forms the **D component** in the insect-based Biodiversity Intactness (IBI) indicator. Predictions are based on insect occurrence records and Earth Observation environmental covariates.

This model supports biodiversity monitoring, conservation planning, land-use impact studies, and policy reporting aligned with the Kunming–Montreal Global Biodiversity Framework.

---

## 🎯 Intended Use

### ✔ Recommended Use Cases
- Mapping insect habitat suitability at continental to national scales  
- Large-scale biodiversity monitoring  
- Conservation prioritization and ecological risk assessments  
- Input layer for the IBI indicator  
- Evaluating land-use and climate impacts on insect diversity  
- Identifying sampling gaps and uncertainty hotspots  

### ✖ Not Designed For
- Predicting species-level presence/absence  
- Population abundance estimation  
- Local-scale modelling (<500 m)  
- Regulatory environmental assessments without additional field data  
- Forecasting under future climate scenarios  

---

## 📊 Model Details

**Model type:** Random Forest (scikit-learn)  
**Output:** Continuous suitability score (0–1)  
**Resolution:** 1 km grid  
**Spatial domain:** Africa (including Madagascar)  
**Version:** v1.2 (2025)

---

## 📥 Training Data

### Occurrence Sources
- Global Biodiversity Information Facility (GBIF)  
- iNaturalist  
- GenBank  

### Taxa Included
- **Coleoptera** (beetles)  
- **Lepidoptera** (butterflies)  
- **Odonata** (dragonflies/damselflies)  

### Cleaning & Filtering
- Removed duplicates  
- Excluded points with >5 km uncertainty  
- Removed ecologically impossible points (urban, permanent water, extreme elevation)  
- Biome-stratified **5 km spatial thinning**  
- Final dataset: **~20,134 high-quality records**  

### Pseudo-Absences
100,000 background points generated in ecologically plausible regions (excluding hyper-arid deserts, permanent water, and urban areas).

---

## 🌍 Predictor Variables

| Variable | Source | Notes |
|---------|--------|-------|
| Brightness (Sentinel-2) | Tasseled Cap | Bare soil, exposure |
| Greenness (Sentinel-2) | Tasseled Cap | Vegetation vigor |
| Wetness (Sentinel-2) | Tasseled Cap | Moisture, canopy composition |
| Canopy Height | GEDI (2020) | 25 m grid, normalized |
| Elevation | SRTM | Gap-filled DEM |

All predictors normalized using z-score transformation.

---

## 🧪 Training Procedure

### Cross-Validation
- **5-fold spatial block CV**
- Blocks ~200 km × 200 km  
- Prevents spatial leakage and inflated performance estimates  

### Final Training
- 70:30 train–test split  
- Hyperparameters tuned with randomized search  

### Hyperparameters
