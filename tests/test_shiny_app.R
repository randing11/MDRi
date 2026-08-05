app <- source(file.path("inst", "shiny", "app.R"), local = TRUE)$value
stopifnot(inherits(app, "shiny.appobj"))
cat("Shiny app construction passed.\n")
