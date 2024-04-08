#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#         Volume calculation            #
#                                       #
#########################################




## Check for the required packages and install them if missing


# The list of required packages
CRAN_packages <- c("tidyverse", "optparse", "this.path", "outliers")


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


# Trim the caliper measurements to all the uCT entries based on date
ctrl_caliper <- list()
temp_df <- data.frame()
for (i in seq_len(length(caliper_measurements))) {
    for (e in seq_len(length(unique(uCT_volumes$Date_calip)))) {
        temp_df <- caliper_measurements[[i]][grep(pattern = unique(uCT_volumes$Date_calip)[e], x = caliper_measurements[[i]]$Dates), ]
        print(temp_df)
        if (nrow(temp_df) > 0) {
        ctrl_caliper[[i]] <- temp_df

    }
    }
    
}


# Trim the caliper measurements to all the uCT entries based on mouse ID
for (i in seq_len(length(unique(uCT_volumes$X)))) {
    if (grepl(pattern = "G1", x = unique(uCT_volumes$X)[i])) {
        ctrl_caliper[[1]] <- ctrl_caliper[[1]][, uCT_volumes$Mouse.ID[grep(pattern = "G1", uCT_volumes$X)]]
    } else {
        ctrl_caliper[[2]] <- ctrl_caliper[[2]][, uCT_volumes$Mouse.ID[grep(pattern = "G2", uCT_volumes$X)]]
    }
}


# Bind the two ctrl caliper DFs together
ctrl_caliper_df <- data.frame()
ctrl_caliper_df <- do.call(cbind, ctrl_caliper)
rownames(ctrl_caliper_df) <- c("L", "W")


# Trim the uCT volumes I will start with
G2_ctrl_uCT <- uCT_volumes[4:7, ]


# Create a unified DF
unified_df <- cbind(uCT_volumes[, c(5, 2, 4)], t(ctrl_caliper_df))




## Calculate the f constants
# NOTE: original formula V = (pi/6) * f * (l * w)^(3/2)
# NOTE: formula solved to f f = V / (pi/6) * (l * w)^(3/2)


# First calculate for a small control set
f_constants <- cbind(G2_ctrl_uCT[, c(5, 2, 4)], t(G2_ctrl_caliper))

for (i in seq_len(nrow(G2_ctrl_uCT))) {
    f_constants$calc_f[i] <- G2_ctrl_uCT$Tumor.volume..mm3.[i] / ((pi / 6) * (G2_ctrl_caliper[1, i] * G2_ctrl_caliper[2, i])^(3 / 2))

}




## Remove the outliers and calculate the mean


# Remove the outliers using Grubb's test (part of the outlier package)
calc_f_values <- f_constants$calc_f

for (i in seq_len(length(calc_f_values))) {
    grubbs_res_lower <- outliers::grubbs.test(calc_f_values)

    if (grubbs_res_lower$p.value < 0.05) {
        lower_outlier <- stringr::str_extract(grubbs_res_lower$alternative, "[0-9]+\\.[0-9]+")
        match_l_indices <- grep(pattern = lower_outlier, x = as.character(calc_f_values))
        
        if (length(match_l_indices) > 0) {
            calc_f_values <- calc_f_values[-match_l_indices]
        }
    } else {
        message("No significant outlier found on the left tail.")

    }

    grubbs_res_upper <- outliers::grubbs.test(calc_f_values, opposite = TRUE)

    if (grubbs_res_upper$p.value < 0.05) {
        upper_outlier <- stringr::str_extract(grubbs_res_upper$alternative, "[0-9]+\\.[0-9]+")
        match_u_indices <- grep(pattern = upper_outlier, x = as.character(calc_f_values))
        
        if (length(match_u_indices) > 0) {
            calc_f_values <- calc_f_values[-match_u_indices]
        }
    } else {
        message("No significant outlier found on the right tail.")

    }

}


# Calculate the mean based on the clean f-values
mean_f <- mean(calc_f_values)




## Re-estimate the tumor volumes with the mean_f for the small ctrl set
for (i in seq_len(nrow(f_constants))) {
    f_constants$est_V[i] <- (pi / 6) * mean_f * (f_constants$L[i] * f_constants$W[i])^(3 / 2)
    #print((pi / 6) * 1 * (f_constants$L[i] * f_constants$W[i])^(3 / 2))

}


# Test the goodness of fit with the estimated Volumes
cont_table <- t(f_constants[, c(3, 7)])
colnames(cont_table) <- f_constants$Mouse.ID

cs <- chisq.test(cont_table)

corr_f <- vector(mode = "numeric")
repeats <- 1
est_V <- vector(mode = "numeric")
if (cs$p.value < 0.05) {
    while (cs$p.value < 0.05) {
        mean_f <- mean_f - 0.1

        for (i in seq_len(nrow(f_constants))) {
            est_V[i] <- (pi/6) * mean_f * (f_constants$L[i] * f_constants$W[i])^(3/2)
            #print((pi / 6) * 1 * (f_constants$L[i] * f_constants$W[i])^(3 / 2))

        }

        cs <- chisq.test(x = rbind(f_constants$Tumor.volume..mm3., est_V))
        print(cs$p.value)
        repeats <- repeats + 1

        if (repeats >= 100) {
            break
        }

    }
}



# Calculate the f for the full used ctrl set
for (i in seq_len(nrow(unified_df))) {
    unified_df$calc_f[i] <- unified_df$Tumor.volume..mm3.[i] / ((pi / 6) * (unified_df$L[i] * unified_df$W[i])^(3 / 2))

}


# Calculate the mean f and re-estimate the tumor volumes with the mean_f for the full ctrl set
full_mean_f <- mean(unified_df$calc_f)

for (i in seq_len(nrow(unified_df))) {
    #f_constants$est_V[i] <- (pi/6) * mean_f * (f_constants$L[i] * f_constants$W[i])^(3/2)
    print((pi / 6) * full_mean_f * (unified_df$L[i] * unified_df$W[i])^(3 / 2))

}










# Calculate the squared distances from the mean
sq_distances <- vector(mode = "numeric", length = nrow(f_constants))
for (i in seq_len(nrow(f_constants))) {
    sq_distances[i] <- (mean_f - f_constants$calc_f[i])^2

}
