cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd, value = TRUE)
root <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  normalizePath(getwd())
}

shiny::runApp(file.path(root, "inst", "shiny"),
              host = "0.0.0.0",
              port = 3839,
              launch.browser = interactive())
