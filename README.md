# FPCA Explorer - Interactive Shiny Application

An interactive R Shiny application for exploring Functional Principal Component Analysis (FPCA/PACE) results. This app provides comprehensive visualization and exploration tools for FPCA bundles saved as `.rds` files.

## Features

### 1. **FPCA Summary Panel**
- Display number of principal components
- Variance table showing eigenvalues, variance proportion, and cumulative variance
- Bar plot of variance proportion across PCs

### 2. **Eigenfunction Visualization**
- Plot eigenfunction shapes on a shared time grid
- Select multiple PCs to display simultaneously
- Customizable Y-axis limits
- Optional mean function overlay

### 3. **PC Effect Plots**
- Visualize how each PC changes curves
- Shows mean ± k·SD along selected PCs
- Adjustable SD range and step size
- Faceted display for multiple PCs
- Mean curve highlighted with thicker line

### 4. **Score Data Table**
- Merged PC scores with metadata
- Searchable and filterable interactive table
- All metadata columns accessible

### 5. **Interactive Score Scatter Plot**
- Choose any two PCs for X and Y axes
- Color and shape aesthetics from metadata
- Interactive tooltips showing key information
- Click points to trigger curve reconstruction
- Built with Plotly for smooth interaction

### 6. **Click-to-Reconstruct Curves**
Two reconstruction modes:

- **Full observation mode**: Reconstruct using first K components (adjustable)
- **Axes-only mode**: Reconstruct using only the two PCs shown on scatter axes

Features:
- Real-time reconstruction on point click
- Optional mean curve overlay
- Clear visualization of reconstruction quality

## Installation

### Required R Packages

```r
install.packages(c(
  "shiny",
  "tidyverse",
  "plotly",
  "DT",
  "funData"
))
```

### Quick Start

1. Clone or download the app file: `fpca_explorer_app.R`

2. Generate sample data (optional, for testing):
```r
source("generate_sample_data.R")
```

3. Run the app:
```r
library(shiny)
runApp("fpca_explorer_app.R")
```

Or in RStudio: Open `fpca_explorer_app.R` and click "Run App"

## FPCA Bundle Specification

Your `.rds` file must contain a named list with the following structure:

### Required Fields

- **`scores`** (numeric matrix, n × npc)  
  Matrix of PC scores, one row per observation

- **`functions`** (funData object, length = npc)  
  Eigenfunctions as a funData object
  - Access eigenfunction matrix: `functions@X` (npc × T)
  - Access time grid: `functions@argvals[[1]]` (length T)

- **`values`** (numeric vector, length = npc)  
  Eigenvalues for each PC

- **`file_order`** (character vector, length = n)  
  Keys/identifiers for each observation (row in scores)

### Optional Fields

- **`fit`** (list or funData)  
  Fitted curves (not required for app functionality)
  - If present, can contain `mu` for mean function

- **`mu`** (funData, 1 × T)  
  Mean function to overlay on plots
  - Can be at top level or inside `fit$mu`

- **`meta`** (data.frame or NULL)  
  Metadata with one row per key
  - Must contain the key column specified in `map_col`
  - Used for coloring/shaping scatter plots

- **`map_col`** (character, default = "File")  
  Name of the key column in `meta` that matches `file_order`

- **`created_at`** (timestamp)  
  Bundle creation timestamp

### Validation

The app validates on load:
- `nrow(scores) == length(file_order)`
- `length(values) == ncol(scores)`
- `functions` contains at least `npc` components
- If `meta` exists: `map_col %in% names(meta)`
- All `file_order` keys exist in `meta[[map_col]]` (warns if missing)

## Example Bundle Creation

```r
library(funData)

# Create eigenfunctions
time_vals <- seq(0, 1, length.out = 200)
phi_mat <- matrix(rnorm(5 * 200), nrow = 5)  # 5 PCs
functions <- funData(argvals = list(time_vals), X = phi_mat)

# Create scores
scores <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)  # 100 obs

# Create bundle
fpca_bundle <- list(
  scores = scores,
  functions = functions,
  values = c(5.0, 3.0, 1.5, 0.8, 0.4),
  file_order = paste0("Subject_", 1:100),
  meta = data.frame(
    File = paste0("Subject_", 1:100),
    Group = sample(c("A", "B"), 100, replace = TRUE),
    Age = rnorm(100, 45, 10)
  ),
  map_col = "File",
  created_at = Sys.time()
)

saveRDS(fpca_bundle, "my_fpca_bundle.rds")
```

## Usage Guide

### Basic Workflow

1. **Upload Data**
   - Click "Upload FPCA Bundle (.rds)" in sidebar
   - Select your `.rds` file
   - App validates and loads the bundle

2. **Explore Summary**
   - View number of PCs and variance explained
   - Check the variance table and bar plot

3. **View Eigenfunctions**
   - Go to "Eigenfunctions" tab
   - Select PCs to display (default: first 3)
   - Adjust Y-limits if needed

4. **Examine PC Effects**
   - Go to "PC Effects" tab
   - Select PCs to visualize
   - Adjust SD range to see curve variations
   - Mean (fraction = 0) shown with thicker line

5. **Explore Scores**
   - Go to "Score Data" tab to see merged table
   - Use filters to search for specific observations

6. **Interactive Analysis**
   - Go to "Interactive Scatter" tab
   - Choose X and Y axis PCs
   - Select color/shape aesthetics from metadata
   - **Click any point** to see reconstruction
   - Toggle between full/axes-only reconstruction
   - Adjust K components for full reconstruction

### Tips

- **Reconstruction Modes**:
  - Use "Full observation" to see how well K components reconstruct the curve
  - Use "Selected axes only" to understand the 2D projection shown in the scatter plot

- **Metadata Aesthetics**:
  - Use color to highlight groups or continuous variables
  - Use shape for categorical variables (works best with <7 categories)

- **Performance**:
  - For large datasets (n > 1000), consider using fewer PCs initially
  - Interactive scatter is optimized with Plotly for smooth interaction

## Customization

The app is designed to be easily customizable. Key sections to modify:

- **Theme**: Change `theme_minimal()` to your preferred ggplot2 theme
- **Colors**: Modify color scales in individual plots
- **Layout**: Adjust `sidebarPanel` width or move controls
- **Validation**: Add custom validation rules in the bundle loading section

## Troubleshooting

### Common Issues

**"Missing required fields" error**
- Ensure your bundle has `scores`, `functions`, `values`, `file_order`

**"Mismatch" errors**
- Check that matrix dimensions align correctly
- Verify `file_order` length matches number of observations

**Plots not updating**
- Check that all required inputs are selected
- Look for error messages in R console

**Click reconstruction not working**
- Make sure you're in the "Interactive Scatter" tab
- Click directly on a point (not empty space)
- Check R console for any error messages

**Mean function not showing**
- Ensure `mu` or `fit$mu` exists in your bundle
- Must be a funData object with matching time grid

## Technical Details

### Dependencies
- **shiny**: App framework
- **tidyverse**: Data manipulation and plotting
- **plotly**: Interactive scatter plot with click events
- **DT**: Interactive data tables
- **funData**: Functional data structures

### Data Flow
1. Bundle loaded and validated reactively
2. Summary computed once on load
3. All plots update reactively based on user inputs
4. Click events trigger reconstruction in real-time
5. No data caching - computation on demand

### Performance Considerations
- Eigenfunction extraction is efficient (direct matrix access)
- Reconstruction computed only for selected point
- Large metadata tables handled by DT pagination
- Plotly downsampling for very large scatter plots (automatic)

## Contact & Support

For issues, questions, or feature requests, please:
- Check the troubleshooting section above
- Review the FPCA bundle specification
- Examine the example data generator code

## License

This app is provided as-is for research and educational purposes.

---

**Version**: 1.0  
**Last Updated**: January 2025  
**Tested with**: R 4.3+, Shiny 1.7+, funData 1.3+
