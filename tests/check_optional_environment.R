cat("R:", R.version.string, "\n")
cat("Library paths:\n")
print(.libPaths())
for (pkg in c("zellkonverter", "basilisk", "reticulate", "UCSCXenaTools")) {
  cat(pkg, requireNamespace(pkg, quietly = TRUE), "\n")
  if (requireNamespace(pkg, quietly = TRUE)) cat("  version:", as.character(packageVersion(pkg)), "\n")
}
if (requireNamespace("zellkonverter", quietly = TRUE)) {
  cat("readH5AD arguments:\n")
  print(args(zellkonverter::readH5AD))
  cat("writeH5AD arguments:\n")
  print(args(zellkonverter::writeH5AD))
}
cat("Environment variables:\n")
print(Sys.getenv(c("R_HOME", "R_LIBS", "R_LIBS_USER", "BASILISK_USE_SYSTEM_DIR", "RETICULATE_PYTHON")))
