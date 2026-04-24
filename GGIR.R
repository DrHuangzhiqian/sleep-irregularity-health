library(GGIR)
library(GGIRread)
args <- commandArgs(trailingOnly = TRUE)
f0 <- as.numeric(args[1])
f1 <- as.numeric(args[2])
# folder_path <- "/public/mig_old_storage/home2/UKB_Preprocessed_Data/Behavior/ukb90001"

cwa_files <- list.files(path = folder_path, pattern = "\\.cwa$", full.names = TRUE, recursive = TRUE)
file_permissions <- sapply(cwa_files, function(x) file.access(x, mode = 4) == 0)
accessible_files <- cwa_files[file_permissions]
GGIR(datadir=accessible_files,outputdir="/public/share/tmp/GGIR_zwc",studyname=paste0(f0, "-", f1),mode=c(1,2,3,4,5),f0=f0,f1=f1)
