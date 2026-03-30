# Tests for SpaCET.deconvolution.bulk
# Mirrors test cases from spatial-gpu/tests/test_deconvolution/test_bulk.py

# Helper: create a minimal SpaCET object with synthetic bulk RNA-seq counts
make_bulk_spacet <- function(n_samples=30, n_genes=200, seed=42) {
  set.seed(seed)

  # Load real reference gene names so gene overlap exists
  ref_path <- system.file("extdata", "combRef_0.5.rda", package = "SpaCET")
  if(!file.exists(ref_path)) return(NULL)
  load(ref_path)
  ref_genes <- rownames(Ref$refProfiles)

  # Use first n_genes reference genes
  if(length(ref_genes) >= n_genes) {
    gene_names <- ref_genes[1:n_genes]
  } else {
    gene_names <- c(ref_genes, paste0("SynGene", seq_len(n_genes - length(ref_genes))))
  }

  # Poisson counts with gene-specific rates
  gene_means <- exp(rnorm(n_genes) * 1.5 + 2)
  counts <- matrix(
    rpois(n_samples * n_genes, rep(gene_means, each=n_samples)),
    nrow=n_genes, ncol=n_samples
  )
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Sample_", formatC(1:n_samples, width=3, flag="0"))

  # Dummy spatial coordinates (required by create.SpaCET.object)
  coords <- data.frame(X = seq_len(n_samples), Y = rep(1, n_samples))
  rownames(coords) <- colnames(counts)

  create.SpaCET.object(
    counts = counts,
    spotCoordinates = coords,
    platform = "OldST"
  )
}


# ---------------------------------------------------------------------------
# Normal tissue branch (mal_prop = 0)
# ---------------------------------------------------------------------------

test_that("bulk deconv runs with cancer_type='normal'", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  expect_s4_class(obj, "SpaCET")
  expect_true(!is.null(obj@results$deconvolution$propMat))
})

test_that("normal tissue sets malProp to zero", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  malProp <- obj@results$deconvolution$malRes$malProp
  expect_true(all(malProp == 0))
})

test_that("normal tissue propMat has correct dimensions", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  pm <- obj@results$deconvolution$propMat
  expect_true(nrow(pm) > 5)  # multiple cell types
  expect_equal(ncol(pm), ncol(obj@input$counts))
})

test_that("normal tissue fractions are non-negative", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  pm <- obj@results$deconvolution$propMat
  expect_true(all(pm >= -1e-10))
})

test_that("normal tissue Malignant row is zero", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  pm <- obj@results$deconvolution$propMat
  expect_true("Malignant" %in% rownames(pm))
  expect_true(all(abs(pm["Malignant",]) < 1e-10))
})


# ---------------------------------------------------------------------------
# External mal_prop branch
# ---------------------------------------------------------------------------

test_that("bulk deconv accepts external malProp array", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  set.seed(99)
  mal <- runif(ncol(obj@input$counts), 0.1, 0.6)
  obj <- SpaCET.deconvolution.bulk(obj, cancerType="BRCA", malProp=mal, coreNo=1)

  expect_true(!is.null(obj@results$deconvolution$propMat))
  expect_true(!is.null(obj@results$deconvolution$malRes$malRef))
})

test_that("bulk deconv accepts named malProp vector", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  set.seed(99)
  mal <- runif(ncol(obj@input$counts), 0.1, 0.6)
  names(mal) <- colnames(obj@input$counts)
  obj <- SpaCET.deconvolution.bulk(obj, cancerType="BRCA", malProp=mal, coreNo=1)

  stored <- obj@results$deconvolution$malRes$malProp
  expect_equal(length(stored), ncol(obj@input$counts))
})

test_that("external malProp clips out-of-range values", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  n <- ncol(obj@input$counts)
  mal <- c(rep(-0.5, n %/% 2), rep(1.5, n - n %/% 2))
  obj <- SpaCET.deconvolution.bulk(obj, cancerType="BRCA", malProp=mal, coreNo=1)

  stored <- obj@results$deconvolution$malRes$malProp
  expect_true(all(stored >= 0))
  expect_true(all(stored <= 1))
})

test_that("higher malProp yields higher Malignant fraction", {
  obj_low <- make_bulk_spacet()
  skip_if(is.null(obj_low), "Reference data not available")
  obj_high <- make_bulk_spacet()

  n <- ncol(obj_low@input$counts)
  obj_low <- SpaCET.deconvolution.bulk(obj_low, cancerType="BRCA", malProp=rep(0.1, n), coreNo=1)
  obj_high <- SpaCET.deconvolution.bulk(obj_high, cancerType="BRCA", malProp=rep(0.8, n), coreNo=1)

  low_mal <- mean(obj_low@results$deconvolution$propMat["Malignant",])
  high_mal <- mean(obj_high@results$deconvolution$propMat["Malignant",])
  expect_true(high_mal > low_mal)
})


# ---------------------------------------------------------------------------
# Signature-based inference branch
# ---------------------------------------------------------------------------

test_that("signature-based bulk deconv runs for BRCA", {
  obj <- make_bulk_spacet(n_samples=30, n_genes=500)
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="BRCA", coreNo=1)

  pm <- obj@results$deconvolution$propMat
  malProp <- obj@results$deconvolution$malRes$malProp

  expect_true(nrow(pm) > 5)
  expect_equal(length(malProp), ncol(obj@input$counts))
  expect_true(all(malProp >= 0))
  expect_true(all(malProp <= 1))
})

test_that("inferMal_bulk returns correct structure", {
  obj <- make_bulk_spacet(n_samples=30, n_genes=500)
  skip_if(is.null(obj), "Reference data not available")

  st <- obj@input$counts
  st <- st[Matrix::rowSums(st) > 0,]

  result <- SpaCET:::inferMal_bulk(st, cancerType="BRCA", signatureType=NULL)

  expect_true(is.list(result))
  expect_true("malRef" %in% names(result))
  expect_true("malProp" %in% names(result))
  expect_equal(length(result$malProp), ncol(st))
  expect_true(all(result$malProp >= 0))
  expect_true(all(result$malProp <= 1))
})

test_that("inferMal_bulk with forced CNA signature", {
  obj <- make_bulk_spacet(n_samples=30, n_genes=500)
  skip_if(is.null(obj), "Reference data not available")

  st <- obj@input$counts
  st <- st[Matrix::rowSums(st) > 0,]

  result <- SpaCET:::inferMal_bulk(st, cancerType="BRCA", signatureType="CNA")

  expect_true(is.list(result))
  expect_equal(length(result$malProp), ncol(st))
})


# ---------------------------------------------------------------------------
# Output structure validation
# ---------------------------------------------------------------------------

test_that("bulk deconv stores results in standard SpaCET slots", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  expect_true(!is.null(obj@results$deconvolution$propMat))
  expect_true(!is.null(obj@results$deconvolution$malRes))
  expect_true(!is.null(obj@results$deconvolution$Ref))
})

test_that("bulk deconv propMat has expected cell types", {
  obj <- make_bulk_spacet()
  skip_if(is.null(obj), "Reference data not available")

  obj <- SpaCET.deconvolution.bulk(obj, cancerType="normal", coreNo=1)

  pm <- obj@results$deconvolution$propMat
  cell_types <- rownames(pm)
  expected <- c("B cell", "T cell", "CAF", "Macrophage", "Endothelial")
  found <- intersect(cell_types, expected)
  expect_true(length(found) >= 2)
})
