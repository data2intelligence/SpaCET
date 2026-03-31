test_that("prepareSpatialCoords returns correct structure for Visium", {
  # Mock visualVector
  visualVector <- c(0.5, 0.3, 0.8)
  names(visualVector) <- c("100x200", "101x201", "102x202")
  spotID <- c("spot1", "spot2", "spot3")

  image <- list(path = NA, grob = NULL)

  result <- SpaCET:::prepareSpatialCoords(
    visualVector, image, "Visium", FALSE, "CaptureArea", NULL, spotID
  )

  expect_true(is.list(result))
  expect_true("fig.df" %in% names(result))
  expect_true("xDiml" %in% names(result))
  expect_true("yDiml" %in% names(result))
  expect_equal(nrow(result$fig.df), 3)
  expect_true("spotID" %in% names(result$fig.df))
})

test_that("prepareSpatialCoords works for non-Visium without image", {
  visualVector <- c(0.5, 0.3)
  names(visualVector) <- c("10x20", "30x40")
  spotID <- c("s1", "s2")

  image <- list(path = NA)

  result <- SpaCET:::prepareSpatialCoords(
    visualVector, image, "OldST", FALSE, "CaptureArea", NULL, spotID
  )

  expect_equal(nrow(result$fig.df), 2)
  expect_null(result$img_raster)
})

test_that("clipLRNetworkScore clips and transforms correctly", {
  mat <- matrix(c(1, 2.0, 0.001, 1, 0.3, 0.1), nrow = 3)

  result <- SpaCET:::clipLRNetworkScore(mat)

  # Row 2 should be clipped to [0.5, 1.5]
  expect_true(all(result[2,] >= 0.5))
  expect_true(all(result[2,] <= 1.5))
  # Row 3 should be -log10 transformed
  expect_equal(result[3, 1], -log10(0.001))
})

test_that("prepareSpatialData returns correct structure for QualityControl", {
  skip("Requires full SpaCET object with QC results")
})

test_that("interactive='plotly' rejects CellTypeComposition", {
  skip("Requires full SpaCET object")
})
