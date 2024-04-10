#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#         Volume calculation            #
#                                       #
#########################################






################################################# MARK: Package and path #################################################
                                                #       management       #




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




## Define the paths for the script and I/O files


# Define the path to the directory this script is in
script_dir_path <- this.path::here()
setwd(script_dir_path)


# Define the path to the intermediate input/output file directory
intermediate_IO <- paste0(script_dir_path, "/", "Intermediate_IO")


################################################# Section end #################################################






################################################# MARK: Data reading #################################################
                                                #   and management   #




## Loading and cleaning the uCT data files needed for the calculations
## NOTE: the first sub-function needs a new file with the uCT measurements 


# This function loads the uCT data from the intermediate I/O folder for further processing
read_uCT_data = function(data_path = intermediate_IO) {
    #Initialize a list onto which the uCT data will be added
    uCT_measurements <- list()
    
    #List the files found in the intermediate I/O folder
    intermediate_IO_files <- list.files(path = data_path)
    
    #Grep the uCT file names for the checking if statement
    uCT_file_names <- intermediate_IO_files[grep(pattern = "uCT", x = intermediate_IO_files)]

    #This main if statement will check if the correct files are in the folder
    if (length(uCT_file_names) > 0) {
        
        #Read each uCT .csv file and add them to the output list
        for (item in intermediate_IO_files) {
            if (grepl(pattern = "uCT", x = item) == TRUE) {
                uCT_measurements[[item]] <- read.csv(file = paste0(data_path, "/", item), sep = ";")
            }
        }

        #Return the uCT data
        return(uCT_measurements)

    } else {
        warning("No uCT input table was found in the ", data_path, " directory.")
    }

}


# This function cleans the uCT data by removing entries without corresponding caliper measurements
clean_uCT_data = function(read_uCT_output_list) {
    #Initialize a new list for the clean uCT dataframes
    clean_uCT_measurement <- list()
    
    #This main for loop will walk through the uCT list of dataframes        
    for (index in seq_len(length(read_uCT_output_list))) {
        #This if statement checks if the column necessary for cleaning is present in the uCT dataframe
        if (is.element(el = "Has_corresponding_caliper_measurements", set = colnames(read_uCT_output_list[[index]])) == TRUE) {
            #Remove entries without corresponding caliper measurements
            clean_uCT_df <- dplyr::filter(.data = read_uCT_output_list[[index]], Has_corresponding_caliper_measurements == "Yes")
            
            #Resets the row numbers
            rownames(clean_uCT_df) <- seq_len(nrow(clean_uCT_df))

            #Add the clean uCT df to the output list
            clean_uCT_measurement[[index]] <- clean_uCT_df

        } else {
            #Print a warning about the missing column
            warning("The column: 'Has_corresponding_caliper_measurements', was not found in the file: ", names(read_uCT_output_list)[index],   
                    "\n  Assuming all entries would be 'Yes' the unmodified uCT_data_frame will be added to the output list")
            
            #Add the unmodified uCT df to the output list
            clean_uCT_measurement[[index]] <- read_uCT_output_list[[index]]

        }
    }    
        
    #Return the clean uCT output list
    return(clean_uCT_measurement)

}




## Loading the processed caliper measurement data files needed for the calculations
## NOTE: this sub-function needs the clean output from the previous module


# This function loads the cleaned output files from the Data-processer module
read_caliper_data = function(data_path = intermediate_IO) {
    #Initialize a list onto which the caliper measurements will be added
    caliper_measurements <- list()

    #List the files found in the intermediate I/O folder
    intermediate_IO_files <- list.files(path = data_path)
    
    #Grep the uCT file names for the checking if statement
    processed_caliper_file_names <- intermediate_IO_files[grep(pattern = "processed", x = intermediate_IO_files)]
    
    #This main if statement will check if the correct files are in the folder
    if (length(processed_caliper_file_names) > 0) {
        #Add the measurements to the list with a for loop
        for (item in intermediate_IO_files) {
            if (grepl("processed", item) == TRUE) {
                caliper_measurements[[item]] <- read.csv(file = paste0(data_path, "/", item))

            }
            
        }
    } else {
        warning("No processed caliper measurement files was found in the ", data_path, " directory.")

    }

    #Name the elements of the caliper_measurements list based on the original processed I/O files
    names(caliper_measurements) <- processed_caliper_file_names

    #Return the caliper_measurement list
    return(caliper_measurements)
    
}


################################################# Section end #################################################






#################################################    MARK: Data processing    #################################################
                                                # for the f-factor estimation #




## Processing the caliper measurements to match the entries of the uCT data


# Trim the caliper measurements to all the uCT entries based on date
# This function will process the caliper measurements to fit the dates and samples recorded in the uCT measurements
fit_caliper_measurements = function(read_caliper_data_output_list, clean_uCT_data_output_list) {
    #Initialize a list onto which the trimmed down caliper measurements will be added
    trimmed_caliper_list <- list()
    
    #Initialize a temporary dataframe on which the processing can take place
    temp_df <- data.frame()

    #initialize a list onto which the unique dates of the corresponding caliper measurements will be added
    corresponding_dates <- vector(mode = "list", length = length(clean_uCT_data_output_list))

    #A for loop which traverses the clean uCT measurement list to assign the corresponding dates to the proper list
    for (i in seq_along(clean_uCT_data_output_list)) {
            
            #An if statements which controls which date columns will be used for sample matching
            if (is.element(el = "Corresponding_caliper_measurement_date", set = colnames(x = clean_uCT_data_output_list[[i]]))) {
                
                #Assign the corresponding dates to the proper list
                corresponding_dates[[i]] <- unique(clean_uCT_data_output_list[[i]]$Corresponding_caliper_measurement_date)
            } else {
                warning("The column: 'Corresponding_caliper_measurement_date' was not found in the", i, "element of the uCT data list.",
                        "\n  Assuming that all caliper measurements correspond with the uCT measurements, the standard Date column will be used.")
                
                #Assign the dates to the proper list
                corresponding_dates[[i]] <- unique(clean_uCT_data_output_list[[i]]$Date)
                
            }
        }




    #The main outer for loop which traverses the caliper measurement list
    for (i in seq_along(read_caliper_data_output_list)) {
        
        #An inner for loop which traverses the corresponding dates list
        for (e in seq_along(corresponding_dates)) {
            
            #The innermost loop which grabs the appropriate dates from the caliper measurements dataframes and assigns them to the temporary dataframe
            for (item in corresponding_dates[[e]]) {
                
                #Assign the date matched caliper measurements to the temp_df
                temp_df <- read_caliper_data_output_list[[i]][grep(pattern = item, x = read_caliper_data_output_list[[i]]$Dates), ]
                    
                    #An if statement to check if rows were assigned
                    if (nrow(temp_df) > 0) {
                        
                        #Assign the temp_df to the trimmed caliper list
                        trimmed_caliper_list[[i]] <- temp_df
                    }
            }
        }
        
    }

    #Name the elements of the caliper_measurement_list base on which treatment group they belong (use the uCT data as reference)



    #Return the trimmed caliper list
    return(trimmed_caliper_list)
    
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


# Create a unified DF
unified_df <- cbind(uCT_volumes[, c(5, 2, 4)], t(ctrl_caliper_df))


# Calculate the f for the full used ctrl set
for (i in seq_len(nrow(unified_df))) {
    unified_df$calc_f[i] <- unified_df$Tumor.volume..mm3.[i] / ((pi / 6) * (unified_df$L[i] * unified_df$W[i])^(3 / 2))

}




## Remove the outliers and calculate the mean


# Remove the outliers using Grubb's test (part of the outlier package)
full_calc_f_values <- unified_df$calc_f

for (i in seq_len(length(full_calc_f_values))) {
    grubbs_res_lower <- outliers::grubbs.test(full_calc_f_values)

    if (grubbs_res_lower$p.value < 0.05) {
        lower_outlier <- stringr::str_extract(grubbs_res_lower$alternative, "[0-9]+\\.[0-9]+")
        match_l_indices <- grep(pattern = lower_outlier, x = as.character(full_calc_f_values))
        
        if (length(match_l_indices) > 0) {
            message("The following lower outliers will be removed: ", match_l_indices)
            full_calc_f_values <- full_calc_f_values[-match_l_indices]
        }
    } else {
        message("No significant outlier found on the left tail.")

    }

    grubbs_res_upper <- outliers::grubbs.test(full_calc_f_values, opposite = TRUE)

    if (grubbs_res_upper$p.value < 0.05) {
        upper_outlier <- stringr::str_extract(grubbs_res_upper$alternative, "[0-9]+\\.[0-9]+")
        match_u_indices <- grep(pattern = upper_outlier, x = as.character(full_calc_f_values))
        
        if (length(match_u_indices) > 0) {
            message("The following lower outliers will be removed: ", match_u_indices)
            full_calc_f_values <- full_calc_f_values[-match_u_indices]
        }
    } else {
        message("No significant outlier found on the right tail.")

    }

}


# Calculate the mean f and re-estimate the tumor volumes with the mean_f for the full ctrl set
full_mean_f <- mean(unified_df$calc_f)

for (i in seq_len(nrow(unified_df))) {
    #f_constants$est_V[i] <- (pi/6) * mean_f * (f_constants$L[i] * f_constants$W[i])^(3/2)
    unified_df$estim_Vol[i] <- (pi / 6) * full_mean_f * (unified_df$L[i] * unified_df$W[i])^(3 / 2)

}

































## Calculate the f constant for the test ctrl uCT measurements


# Trim the caliper measurements to the entries I will start with
G2_ctrl_caliper <- caliper_measurements[[2]][grep(pattern = "21.02.2024", x = caliper_measurements[[2]]$Dates), uCT_volumes$Mouse.ID[4:7]]
rownames(G2_ctrl_caliper) <- c("L", "W")

# Trim the uCT volumes I will start with
G2_ctrl_uCT <- uCT_volumes[4:7, ]




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




## Test the whole process with all of the late control measurements


# Add the actual caliper measurement dates to the uCT_volumes DF
uCT_volumes <- mutate(.data = uCT_volumes,
                Date_calip = c("03.03.2024", "03.03.2024", "03.03.2024", "21.02.2024", "21.02.2024", "21.02.2024", "21.02.2024"),
                .before = Comments)

# Trim the caliper measurements to all the uCT entries based on date
ctrl_caliper <- list()
temp_df <- data.frame()
for (i in seq_len(length(caliper_measurements))) {
    for (e in seq_len(length(unique(uCT_volumes$Date_calip)))) {
        temp_df <- caliper_measurements[[i]][grep(pattern = unique(uCT_volumes$Date_calip)[e], x = caliper_measurements[[i]]$Dates), ]
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


# Create a unified DF
unified_df <- cbind(uCT_volumes[, c(5, 2, 4)], t(ctrl_caliper_df))


# Calculate the f for the full used ctrl set
for (i in seq_len(nrow(unified_df))) {
    unified_df$calc_f[i] <- unified_df$Tumor.volume..mm3.[i] / ((pi / 6) * (unified_df$L[i] * unified_df$W[i])^(3 / 2))

}




## Remove the outliers and calculate the mean


# Remove the outliers using Grubb's test (part of the outlier package)
full_calc_f_values <- unified_df$calc_f

for (i in seq_len(length(full_calc_f_values))) {
    grubbs_res_lower <- outliers::grubbs.test(full_calc_f_values)

    if (grubbs_res_lower$p.value < 0.05) {
        lower_outlier <- stringr::str_extract(grubbs_res_lower$alternative, "[0-9]+\\.[0-9]+")
        match_l_indices <- grep(pattern = lower_outlier, x = as.character(full_calc_f_values))
        
        if (length(match_l_indices) > 0) {
            message("The following lower outliers will be removed: ", match_l_indices)
            full_calc_f_values <- full_calc_f_values[-match_l_indices]
        }
    } else {
        message("No significant outlier found on the left tail.")

    }

    grubbs_res_upper <- outliers::grubbs.test(full_calc_f_values, opposite = TRUE)

    if (grubbs_res_upper$p.value < 0.05) {
        upper_outlier <- stringr::str_extract(grubbs_res_upper$alternative, "[0-9]+\\.[0-9]+")
        match_u_indices <- grep(pattern = upper_outlier, x = as.character(full_calc_f_values))
        
        if (length(match_u_indices) > 0) {
            message("The following lower outliers will be removed: ", match_u_indices)
            full_calc_f_values <- full_calc_f_values[-match_u_indices]
        }
    } else {
        message("No significant outlier found on the right tail.")

    }

}


# Calculate the mean f and re-estimate the tumor volumes with the mean_f for the full ctrl set
full_mean_f <- mean(unified_df$calc_f)

for (i in seq_len(nrow(unified_df))) {
    #f_constants$est_V[i] <- (pi/6) * mean_f * (f_constants$L[i] * f_constants$W[i])^(3/2)
    unified_df$estim_Vol[i] <- (pi / 6) * full_mean_f * (unified_df$L[i] * unified_df$W[i])^(3 / 2)

}













# Calculate the squared distances from the mean
sq_distances <- vector(mode = "numeric", length = nrow(f_constants))
for (i in seq_len(nrow(f_constants))) {
    sq_distances[i] <- (mean_f - f_constants$calc_f[i])^2

}
