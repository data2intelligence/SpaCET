test_that("create.SpaCET.object.10X loads Visium BC data", {
  visium_path <- file.path(system.file("extdata", package = "SpaCET"), "Visium_BC")
  skip_if(!dir.exists(visium_path), "Visium_BC data not available")

  obj <- create.SpaCET.object.10X(visiumPath = visium_path)

  expect_s4_class(obj, "SpaCET")
  expect_true(nrow(obj@input$counts) > 0)
  expect_true(ncol(obj@input$counts) > 0)
  expect_true("spotCoordinates" %in% names(obj@input))
  expect_equal(obj@input$organism, "human")
})

test_that("create.SpaCET.object.10X accepts organism parameter", {
  visium_path <- file.path(system.file("extdata", package = "SpaCET"), "Visium_BC")
  skip_if(!dir.exists(visium_path), "Visium_BC data not available")

  obj <- create.SpaCET.object.10X(visiumPath = visium_path, organism = "mouse")
  expect_equal(obj@input$organism, "mouse")
})

test_that("create.SpaCET.object validates input dimensions", {
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("Gene", 1:10)
  colnames(counts) <- paste0("Spot", 1:10)

  coords <- data.frame(X = runif(10), Y = runif(10))
  rownames(coords) <- paste0("Spot", 1:10)

  obj <- create.SpaCET.object(
    counts = counts,
    spotCoordinates = coords,
    platform = "OldST"
  )

  expect_s4_class(obj, "SpaCET")
  expect_equal(ncol(obj@input$counts), 10)
  expect_equal(nrow(obj@input$counts), 10)
})

test_that("create.SpaCET.object rejects mismatched dimensions", {
  counts <- matrix(1, nrow = 5, ncol = 10)
  colnames(counts) <- paste0("Spot", 1:10)
  rownames(counts) <- paste0("Gene", 1:5)

  coords <- data.frame(X = runif(8), Y = runif(8))
  rownames(coords) <- paste0("Spot", 1:8)

  expect_error(
    create.SpaCET.object(counts = counts, spotCoordinates = coords, platform = "OldST"),
    "identical"
  )
})

test_that("SpaCET.quality.control filters spots", {
  visium_path <- file.path(system.file("extdata", package = "SpaCET"), "Visium_BC")
  skip_if(!dir.exists(visium_path), "Visium_BC data not available")

  obj <- create.SpaCET.object.10X(visiumPath = visium_path)
  obj <- SpaCET.quality.control(obj, min.genes = 1)

  expect_true(!is.null(obj@results$metrics))
  expect_equal(nrow(obj@results$metrics), 2)  # UMI and Gene rows
})

test_that("addVisiumMicrometerCoords adds coordinate columns", {
  coords <- data.frame(
    array_row = c(0, 1, 2),
    array_col = c(0, 1, 2),
    row.names = c("0x0", "1x1", "2x2")
  )

  result <- SpaCET:::addVisiumMicrometerCoords(coords)

  expect_true("coordinate_x_um" %in% names(result))
  expect_true("coordinate_y_um" %in% names(result))
  expect_equal(result[1, "coordinate_x_um"], 0)
})

test_that("mouse2human_mat converts mouse genes", {
  skip_if(!file.exists(system.file("extdata", "Mouse2Human_filter.csv", package = "SpaCET")),
          "Mouse mapping file not available")

  # Create small matrix with known mouse genes
  m2h <- read.csv(system.file("extdata", "Mouse2Human_filter.csv", package = "SpaCET"), row.names = 1)
  mouse_genes <- head(m2h$mouse, 5)

  mat <- matrix(1:10, nrow = 5, ncol = 2)
  rownames(mat) <- mouse_genes

  result <- SpaCET:::mouse2human_mat(mat)

  expect_true(nrow(result) > 0)
  expect_equal(ncol(result), 2)
  # All row names should be human genes
  expect_true(all(rownames(result) %in% m2h$human))
})
