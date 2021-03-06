dt <- list(); dt$R <- paste0(Sys.getenv("OneDriveCommercial"), "/FE_Methoden/", "Allgemein/R_dt_project/")
source(paste0(dt$R,"R/source_spc_files.R"))
source(paste0(dt$R,"R/source_pls.R"))
source(paste0(dt$R,"R/source_read.R"))

# general parameters ####
dt$para$customer = "PepsiCo"
dt$para$beverage = "Pepsi_Light"

setwd(paste0(dt$wd <- paste0(wd$fe$Pepsi$Mastermodelle, dt$para$beverage)))
setwd( print( this.path::this.dir() ) )
# setwd("..")
dt$wd.git <- print( getwd() )

dt$para$location = "Nieder_Roden"
dt$para$line = "L3_PET_CSD"
dt$para$main = paste0(dt$para$beverage, " in ", dt$para$location, ", line ", dt$para$line)
dt$para$model.date <- c("220512")
dt$para$model.pl <- c("00300")
dt$para$wl1 <- c(190)
dt$para$wl2 <- c(598)
dt$para$wl[[1]] <- seq(dt$para$wl1, dt$para$wl2, 1)

# rename R files (run only once)
dt$para$Rfiles <- list.files(getwd(), pattern = ".R$", recursive = F)
file.rename(dt$para$Rfiles, gsub("beverage", dt$para$beverage, dt$para$Rfiles))

