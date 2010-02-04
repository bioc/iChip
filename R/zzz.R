.First.lib <- function(lib, pkg) {
   library.dynam("iChip", pkg, lib)
   cat("iChip 1.4.0 loaded\n")
}

