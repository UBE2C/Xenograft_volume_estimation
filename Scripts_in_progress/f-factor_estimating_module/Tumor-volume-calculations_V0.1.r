#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#         Volume calculation            #
#                                       #
#########################################




## Check for the required packages and install them if missing


# The list of required packages
CRAN_packages <- c("tidyverse", "optparse", "this.path")


# A function to check and install packages
package_controller = function(packages = CRAN_packages) {
    
    #The main if statement to make sure there are packages to check
    if (!is.null(packages)) {
        
        #This loop goes through the package list and feeds them to the checking if statement
        for (element in packages) {
            
            #This if statement checks if the given package is part of the installed packages
            if (!is.element(element, installed.packages())) {
                
                message(paste0("The following package is not installed: ", element, "\n", 
                "Installing ", element, "now"))
                install.packages(element)

            } else {
                message(paste0("The requested package: ", element, "is already installed."))
            }
        }
    } else {
        message("No package was requested.")
    }
    
}


# A function to load required packages
package_loader = function(packages = CRAN_packages) {
    lapply(packages, library, character.only = TRUE)
}




## Define the path to the directory this script is in
script_dir_path <- this.path::here()
setwd(script_dir_path)




## Define the path to the intermediate input file directory
interm_inputs <- paste0(script_dir_path, "/", "Intermediate_input_files")




## Loading and cleaning the data files needed for the calculations


# Loading the control uCT measurements
uCT_volumes <- read.csv(file = paste0(interm_inputs, "/", "uCT_ctrl_measurements.csv"),
                        sep = ";")

# Remove entries without caliper measurements
uCT_volumes <- uCT_volumes[-c(1, 3, 9), ]
rownames(uCT_volumes) <- seq_len(nrow(uCT_volumes))


# Loading the cleaned output files from the Data-processer 
caliper_measurements <- list()
for (element in list.files(interm_inputs)) {
    if (grepl("processed", element) == TRUE) {
        caliper_measurements[[element]] <- read.csv(file = paste0(interm_inputs, "/", element))

    }
    
}




## Calculate the f constant for the test ctrl uCT measurements


# Trim the caliper measurements to the entries I will start with
G2_ctrl_caliper <- caliper_measurements[[2]][grep(pattern = "21.02.2024", x = caliper_measurements[[2]]$Dates), uCT_volumes$Mouse.ID[4:7]]
rownames(G2_ctrl_caliper) <- c("L", "W")


# Trim the uCT volumes I will start with
G2_ctrl_uCT <- uCT_volumes[4:7, ]


# Calculate the f constants
# NOTE: original formula V = (pi/6) * f * (l * w)^(3/2)
# NOTE: formula solved to f f = V / (pi/6) * (l * w)^(3/2)
f_constants <- cbind(G2_ctrl_uCT[, c(5, 2, 4)], t(G2_ctrl_caliper[, c(1, 2)]))

for (i in seq_len(nrow(G2_ctrl_uCT))) {
    f_constants$calc_f[i] <- G2_ctrl_uCT$Tumor.volume..mm3.[i] / ((pi/6) * (G2_ctrl_caliper[1, i] * G2_ctrl_caliper[2, i])^(3/2))

}


# Calculate the mean f and re-estimate the tumor volumes with the mean_f
