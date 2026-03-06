# FPCA Visualiser

An interactive R Shiny application for exploring Functional Principal Component Analysis (FPCA) results.

The app can be accessed online: https://k-li.shinyapps.io/InterVisual-fpca/

## Features

- **Summary**: eigenvalues, variance proportions, and cumulative variance table and bar plot
- **Eigenfunctions**: plot eigenfunction shapes with optional mean curve overlay
- **PC Effect Plots**: visualise how each PC shifts the mean curve (mean ± k·SD)
- **Score Scatter**: interactive scatter plot of any two PCs, coloured/shaped by metadata; click points to select observations
- **Reconstruction**: reconstruct selected curves using either the K leading PCs or only the two displayed axes

## Usage

Upload an FPCA bundle (`.rds`) via the sidebar, or click **Load Example Dataset** to explore with the built-in example. See the **Instruction** tab inside the app for details on preparing your own bundle.

---

**Tested with**: R 4.3+, Shiny 1.7+, funData 1.3+
