dt <- list(); dt$R <- paste0(Sys.getenv("OneDriveCommercial"), "/FE_Methoden/", "Allgemein/R_dt_project/")
source(paste0(dt$R,"R/source_spc_files.R"))
source(paste0(dt$R,"R/source_pls.R"))
source(paste0(dt$R,"R/source_read.R"))

# general parameters ####
dt$para$customer = "PepsiCo"
dt$para$beverage = "Pepsi_Light"

setwd(paste0(dt$wd <- paste0(wd$fe$Pepsi$Mastermodelle, dt$para$beverage)))
setwd( print( this.path::this.dir() ) )
setwd("..")
dt$wd.git <- print( getwd() )

dt$para$location = "Nieder_Roden"
dt$para$line = "L3_PET_CSD"
dt$para$main = paste0(dt$para$beverage, " in ", dt$para$location, ", line ", dt$para$line)
dt$para$model.date <- c("220318")
dt$para$model.pl <- c("00300")
dt$para$wl1 <- c(190)
dt$para$wl2 <- c(598)
dt$para$wl[[1]] <- seq(dt$para$wl1, dt$para$wl2, 1)

dt$para$substance <- c("TA",	"Coffein",	"Aspartam",	"Acesulfam", "Benzoat")
dt$para$unit <- c( bquote("mg L"-1),  bquote("mg L"-1),  bquote("mg L"-1),  bquote("mg L"-1) ,  bquote("mg L"-1))
dt$para$ylab <- c(bquote("TA in mg / 100mL"^-1), bquote("Coffein in mg / L"^-1), bquote("Aspartam in mg / L"^-1), bquote("Acesulfam in mg / L"^-1), bquote("Benzoat in mg / L"^-1))
dt$para$mop.date <- "220522"
dt$para$SOLL <- c(12.84, 120, 340.09, 98.65, 143.92)

# Model Matrix ####
setwd(dt$wd)
setwd("./Modellerstellung")
setwd(paste0("./", dt$para$model.date[1], "_", dt$para$model.pl[1]))

dt$model <- read.csv2( print(grep( "match.csv", dir(), value = T)), dec = ",", sep = ";")
# colnames(dt$model)[colnames(dt$model) %in% "SK"] <- "TA"
# dt$model$Probe <- gsub("SK", "TA", dt$model$Probe)
# dt$model$Probe_Anteil <- gsub("SK", "TA", dt$model$Probe_Anteil)

for(i in 1:length(dt$para$substance)){
  dt$model[ , colnames(dt$model) %in% dt$para$substance[i]] <- dt$model[ , colnames(dt$model) %in% dt$para$substance[i]] * dt$para$SOLL[i] / 100
}
# 

dt$model <- dt$model[ , - which( as.numeric(gsub("X", "", colnames(dt$model))) < 190 | as.numeric(gsub("X", "", colnames(dt$model))) > 598)]

dt$SL <- dt$model[which(dt$model$Probe_Anteil == "SL") , ]
dt$model <- dt$model[which(dt$model$Probe_Anteil != "SL") , ]

# Linearity
setwd(dt$wd)
setwd("./Modellvalidierung")
setwd("./Linearitaet")
dt$lin$raw <- read.csv2( "220524_Linearitaet_Pepsi_Light_Aspartam_Acesulfam_Coffein_Benzoat_Acid.csv" , sep = "\t")
dt$lin$raw <- dt$lin$raw[ order(dt$lin$raw$SOLL) , ]
dt$lin$trs <- .transfer_csv(dt$lin$raw)

# rename R files (run only once)
# dt$para$Rfiles <- list.files(getwd(), pattern = ".R$", recursive = T)
# file.rename(dt$para$Rfiles, gsub("beverage", dt$para$beverage, dt$para$Rfiles))

