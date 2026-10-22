# Installs missing dependencies into the project, without altering the global library.
arg <- grep("^--file=",commandArgs(FALSE),value=TRUE)
root <- if(length(arg)) dirname(dirname(normalizePath(sub("^--file=","",arg[1])))) else normalizePath(".")
lib <- file.path(root,".R-library"); dir.create(lib,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
packages <- c("data.table","jsonlite","digest","urca")
needed <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(needed)) install.packages(needed,lib=lib,repos="https://cloud.r-project.org")
stopifnot(all(vapply(packages,requireNamespace,logical(1),quietly=TRUE)))
cat("R dependencies ready.\n")
