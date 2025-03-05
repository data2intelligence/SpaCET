## SpaCET 1.3.0

*2025-03-05*

* Add cell references of hepatocytes and cholangiocytes so that SpaCET can infer their fractions when deconvolving liver cancer samples.
* Update `SpaCET.deconvolution` to deconvolve normal tissues by setting`adjacentNormal = TRUE` to skip the malignant cell prediction.
* Update `SpaCET.visualize.spatialFeature` to visualize secreted protein signaling activity and pattern.

## SpaCET 1.2.0

*2024-07-06*

* Add a new function `SpaCET.GeneSetScore` to calculate gene set score.
* Update `convert.Seurat` and `addTo.Seurat` to be compatible with Seurat v5.

## SpaCET 1.1.0

*2023-01-30*

* Add an interactive visualization panel to browse the deconvolution and analysis results `SpaCET.visualize.spatialFeature(SpaCET_obj,interactive=TRUE)`.
* Add a new function to combine both the tumor-immune interface and interaction spots `SpaCET.combine.interface()`.

## SpaCET 1.0.0

*2022-12-20*

* Release SpaCET package.
