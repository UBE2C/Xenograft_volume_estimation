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
CRAN_packages <- c("tidyverse", "optparse", "this.path", "outliers", "ggpubr")


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

        #Name the return list elements based on the original file names
        names(uCT_measurements) <- uCT_file_names

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

    #Name the elements of the return list based on the elements of the input list
    names(clean_uCT_measurement) <- names(read_uCT_output_list)
        
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
    ##Define the variables used by this function
    
    #Initialize a list onto which the trimmed down caliper measurements will be added
    trimmed_caliper_list <- list()
    
    #Initialize a temporary dataframe on which the processing can take place
    temp_df <- data.frame()

    #Initialize a list onto which the unique dates of the corresponding caliper measurements will be added
    corresponding_dates <- vector(mode = "list", length = length(clean_uCT_data_output_list))

    #Define the treatment group_IDs shared by both uCT and caliper_measurement files
    group_IDs <- vector(mode = "list", length = length(clean_uCT_data_output_list))

    for (i in seq_along(clean_uCT_data_output_list)) {
        group_IDs[[i]] <- unique(clean_uCT_data_output_list[[i]]$Treatment_group_ID)
    }
    

    ##Start the processing steps

    #First, trim the caliper measurements by date to match the uCT measurements
    
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
    
    #Assign the corresponding_dates to the parent main() environment
    assign("corresponding_dates_main", corresponding_dates, envir = parent.frame())

    #The main outer for loop which traverses the corresponding dates list
    for (i in seq_along(corresponding_dates)) {
        
        #An inner for loop which traverses the elements inside corresponding dates list
        for (e in seq_along(corresponding_dates[[i]])) {
            
            #The innermost loop which traverse the caliper measurements list 
            for (df_i in seq_along(read_caliper_data_output_list)) {
                
                #Assign the date matched caliper measurements to the temp_df
                temp_df <- read_caliper_data_output_list[[df_i]][grep(pattern = corresponding_dates[[i]][[e]], x = read_caliper_data_output_list[[df_i]]$Dates), ]
                
                    #An if statement to check if rows were assigned
                    if (nrow(temp_df) > 0) {
                        
                        #Append the trimmed caliper list with the temp_df 
                        trimmed_caliper_list <- append(x = trimmed_caliper_list, values = list(temp_df))
                        
                    }
            }
        }
        
    }

    #Second, trim the caliper measurements by experimental group and subject to match the uCT measurements
    #Further subset the trimmed_caliper_list elements based on the sample names of the uCT measurements
    #NOTE: this is a complicated loop system so I will add extra annotation to facilitate better understanding

    #The outermost loop which traverses the uCT data based on which we are subsetting, this will provide the uct and the group_ID list index (as they are supposed to be the same)
    for (i in seq_along(clean_uCT_data_output_list)) {
        
        #The first inner loop which traverses the group_IDs list. The unique shared IDs should be stored here for the grep function
        #this will allow to select the right list element(i) and the right identifier for that uCT ID vector (e)
        for (e in seq_along(group_IDs[[i]])) {
            
            #The innermost loop which traverses the trimmed_caliper_measurements list. The Treatment_group_ID column of each list element should contain the 
            #shared IDs so the grep function can identify in which list element the given ID can be found.
            for (df_e in seq_along(trimmed_caliper_list)) {
                
                #An if statement ensuring that if there are multiple list elements to process, and the first element was processed and the column
                #Treatment_group_ID is missing because of that, the whole function won't break but jumps to the next list element
                if (is.element(el = "Treatment_group_ID", colnames(trimmed_caliper_list[[df_e]]))) {

                    #An if statement ensuring that if the given ID is not found in the Treatment_group_ID of a caliper list element
                    #the loop does not break but jumps to the next element
                    if (grepl(pattern = group_IDs[[i]][[e]], x = unique(trimmed_caliper_list[[df_e]]$Treatment_group_ID))) {
                        
                        #The actual part responsible for the dataframe splitting.Consist of two main pieces -
                        #- Piece 1-
                        #clean_uCT_data_output_list[[i]]$Mouse_ID[grep(pattern = group_IDs[[i]][[e]], x = clean_uCT_data_output_list[[i]]$Treatment_group_ID)]
                        #This piece will look for the shared ID in the uCT data shared ID column, and uses the returned row indexes to subset the uCT dataframe's
                        #mouse ID column into the matching mouse IDs
                        #- Piece 2-
                        #trimmed_caliper_list[[df_i]][, -Piece 1-]
                        #This piece will simply split the trimmed caliper measurement dataframes down to the selected columns, returned by -Piece 1-
                        trimmed_caliper_list[[df_e]] <- trimmed_caliper_list[[df_e]][, clean_uCT_data_output_list[[i]]$Mouse_ID[grep(pattern = group_IDs[[i]][[e]], x = clean_uCT_data_output_list[[i]]$Treatment_group_ID)]]
                    } else {
                        next
                    }
                } else {
                    next
                }

                
            }
            
        }
        
    }

    #Return the trimmed caliper list
    return(trimmed_caliper_list)
    
}

#Bind the trimmed caliper measurement dataframes together by columns (animal IDs) and create a clean unified df using the uCT measurement dfs
#This function will carry out the column bind and organizes the dataframes according to the uCT measurements list, then creates a list of unified dataframes
#using the selected column of the cuCT measurements and the transposed bound caliper measurements
bind_and_unify_measurements = function(fit_caliper_measurements_output_list, clean_uCT_data_output_list) {
    
    
    ##Define the variables used in the function

    #Assign the trimmed_caliper_measurement list to a variable which will be trimmed progressively
    shrinking_list <- list()
    shrinking_list <- fit_caliper_measurements_output_list

    #Initialize a new list which will contain the column bound dataframes. Each element will correspond to the
    #elements of the uCT measurements list
    bound_df_list <- list()

    #Initialize an empty vector to store the number of dates in each corresponding date element
    no_of_dates <- vector(mode = "numeric")

    #Initialize a unified dataframe list, which will unify the bound caliper measurements and selected columns from the clean_uCT_list
    unified_df_list <- list()


    ##Start the processing steps

    #This for loop will traverse the corresponding dates list. The length of each element will correspond to the number of dataframes 
    #the trimmed_caliper_measurements list has, so they can be grouped based on the uCT files
    for (i in seq_along(corresponding_dates_main)) {
        
        #Assign the length of the corresponding date list element
        no_of_dates <- length(corresponding_dates_main[[i]])
        
        #Assign the column bound dataframes to the bound_df_list
        bound_df_list[[i]] <- cbind(shrinking_list[[1]], shrinking_list[[no_of_dates]])
        
        #Shrink the caliper measurement df list by the elements bound beforehand
        shrinking_list <- shrinking_list[- c(1, no_of_dates)]
    
    }

    #Set the rownames of the new dataframes
    for (i in seq_along(bound_df_list)){
        rownames(bound_df_list[[i]]) <- c("L", "W")
    }

    #Unify the bound caliper measurements and selected columns from the clean_uCT_list
    for (i in seq_along(clean_uCT_data_output_list)) {
        unified_df_list[[i]] <- cbind(clean_uCT_data_output_list[[i]][, c("Date", "Mouse_ID", "Tumor_volume_.mm3.")], t(bound_df_list[[i]]))
    }


    #Return the bound_df_list
    return(unified_df_list)

}


################################################# Section end #################################################





################################################    MARK: f-constant    #################################################
                                                #  calculation and eval #




## Calculate the f-constants to each uCT-measurement-caliper measurement group

# This function calculates the f-constant for the used uCT measurement set
calculate_f_constants = function(bind_and_unify_measurements_output_list) {
    

    ##Define the variables used in the function

    #Assign the bind_and_unify_measurements_output_list to a variable which will be returned
    measurements_with_f_constants_list <- bind_and_unify_measurements_output_list

    #This outer loop traverses the unified measurement dfs list
    for (i in seq_along(bind_and_unify_measurements_output_list)) {
        
        #This inner loop traverses the unified dataframes themselves and adds a new column "f_constant" and fills it with the
        #calculated f-constants to each sample
        for (e in seq_len(nrow(bind_and_unify_measurements_output_list[[i]]))) {
            measurements_with_f_constants_list[[i]]$f_constants[e] <- bind_and_unify_measurements_output_list[[i]]$Tumor_volume_.mm3.[e] / ((pi / 6) * (bind_and_unify_measurements_output_list[[i]]$L[e] * bind_and_unify_measurements_output_list[[i]]$W[e])^(3 / 2))

        }
    }

    #Return the output list containing the f-constants
    return(measurements_with_f_constants_list)
}


################################################# Section end #################################################





################################################      MARK: Outlier calculation      #################################################
                                                #  and removal among the f-constants #
# MARK: CONTINUE FROM HERE



## Determine if the data is normally distributed
is_data_normal = function(calculate_f_constants_output_list) {
    

    ## Declare dynamic variables
    shapiro_results <- vector(mode = "list", length = length(calculate_f_constants_output_list))


    ## Do the normality test

    # This for loop will carry out a Shapiro-Wilk normality test
    for (i in seq_along(calculate_f_constants_output_list)) {
        shapiro_results[[i]] <- shapiro.test(calculate_f_constants_output_list[[i]]$f_constants)
    }

    # Return the shapiro test results list
    return(shapiro_results)

}




## Non-parametric outlier test - Numeric outlier test (this seems to be the more safe alternative as it is a standard test)
## The function identifies adn removes outliers from the input df based on the f-constants
remove_outlier_f_const_nonparam = function(non_normal_list_element) {
    ## Declare a list to modify
    list_element_to_modify <- non_normal_list_element

    # Declare the output list
    f_constant_outlier_free_measurements <- data.frame(matrix(nrow = nrow(list_element_to_modify), ncol = ncol(list_element_to_modify)))

    ## Declare variables

    # Declare the quartile variables
    Q1 <- vector(mode = "numeric", length = 1)
    Q3 <- vector(mode = "numeric", length = 1)

    # Declare a vector for the IQR - inter quartile range
    IQR <- vector(mode = "numeric", length = 1)

    # Declare vectors containing the upper and lower boundaries
    upper_boundary <- vector(mode = "numeric", length = 1)
    lower_boundary <- vector(mode = "numeric", length = 1)

    # Declare vectors for the upper and lower outliers
    upper_outlier <- vector(mode = "numeric")
    lower_outlier <- vector(mode = "numeric")

    # Order the constants in am ascending order
    ordered_constants <- list_element_to_modify$f_constants[order(list_element_to_modify$f_constants)]


    # Calculate the quartile ranges and allocate them to the appropriate list elements
    Q1 <- quantile(ordered_constants, 0.25, na.rm = TRUE)
    Q3 <- quantile(ordered_constants, 0.75, na.rm = TRUE)

    # Calculate the IQR 
    IQR <- Q3 - Q1

    # Calculate the upper and lower boundaries
    upper_boundary <- Q3 + (1.5 * IQR)
    lower_boundary <- Q1 - (1.5 * IQR)

    if (lower_boundary < 0) {
        lower_boundary <- 0
    }

    # Identify the upper outliers
    upper_outlier <- ordered_constants[ordered_constants > upper_boundary]

    # Identify the lower outliers
    lower_outlier <- ordered_constants[ordered_constants < lower_boundary]
        
    # An if statement to check if there were any upper outliers and if yes to remove them
    if (length(upper_outlier) == 0) {
        message("No outlier found on the right (upper) tail for input list element.")
    } else {
        message("The following elements are outliers on the right (upper) tail, and will be removed: \n", upper_outlier)
            
        #Remove the upper outliers
        list_element_to_modify <- list_element_to_modify[!list_element_to_modify$f_constants %in% upper_outlier, ]
    }

    # An if statement to check if there were any lower outliers and if yes to remove them
    if (length(lower_outlier) == 0) {
        message("No outlier found on the left (lower) tail for input list element.")
    } else {
        message("The following elements are outliers on the left (lower) tail, and will be removed: \n", lower_outlier)

        #Remove the lower outliers
        list_element_to_modify <- list_element_to_modify[!list_element_to_modify$f_constants %in% lower_outlier, ]
    }



    # Assign the f-constant outlier free list element to the output list
    f_constant_outlier_free_measurements <- list_element_to_modify

    # Return the output list
    return(f_constant_outlier_free_measurements)
    
}




## Modified z-score test (non-parametric) - this test is basically the Grubb's test, just with the median and MAD instead of mean and SD
## The ida was taken from the following publication - DOI: (https://doi.org/10.1515/dema-2021-0041), however the publication used the SD with the median
## and not the MAD
## The function identifies adn removes outliers from the input df based on the f-constants
modified_Zscore_test = function(non_normal_list_element, left_tail = TRUE) {
    
    
    ## Define function variables

    # Define an output dataframe
    output_df <- data.frame(matrix(nrow = nrow(non_normal_list_element), ncol = ncol(non_normal_list_element)))

    # Define the output list
    output_list <- vector(mode = "list", length = 2)
    
    # Define the f_constant vector for sorting
    ordered_f_constants <- vector(mode = "numeric", length = nrow(non_normal_list_element))
    
    # Define the median f-constant
    f_const_median <- vector(mode = "numeric", length = 1)

    # Define the median absolute deviation (MAD)
    f_const_mad <- vector(mode = "numeric", length = 1)

    # Define a vector for the z-scores
    z_scores <- vector(mode = "numeric", length = nrow(non_normal_list_element))

    # Define a vector for the minimum z-score
    min_z_score <- vector(mode = "numeric", length = 1)

    # Define a vector for the maximum z-score
    max_z_score <- vector(mode = "numeric", length = 1)

    # Define the critical z-value
    crit_z <- vector(mode = "numeric", length = 1)

    # Max outlier identity
    max_outlier_ident <- vector(mode = "numeric", length = 1)

    # Max outlier index
    max_outlier_index <- vector(mode = "numeric", length = 1)

    # Min outlier identity
    min_outlier_ident <- vector(mode = "numeric", length = 1)

    # Min outlier index
    min_outlier_index <- vector(mode = "numeric", length = 1)

    # Variable to track if no outliers are found
    while_loop_run <- TRUE


    ## Calculate the defined variables

    # Assign and sort the f_constants
    ordered_f_constants <- non_normal_list_element$f_constants[order(non_normal_list_element$f_constants)]

    #Assign the input dataframe as an output dataframe so in case of no outliers the original df will be returned
    output_df <- non_normal_list_element

    # Calculate the median of the f-constants
    f_const_median <- median(ordered_f_constants, na.rm = TRUE)

    # Calculate the MAD
    f_const_mad <- mad(x = non_normal_list_element$f_constants, na.rm = TRUE)

    # Calculate the minimum z-score
    min_z_score <- (f_const_median - min(non_normal_list_element$f_constants)) / f_const_mad

    # Calculate all the z-scores
    for (i in seq_along(non_normal_list_element$f_constants)) {
        z_scores[i] <- abs((non_normal_list_element$f_constants[i] - f_const_median) / f_const_mad)
    
    }

    # Calculate the maximum z-score
    max_z_score <- (max(non_normal_list_element$f_constants) - f_const_median) / f_const_mad

    # Calculate the critical z-score
    crit_z <- (nrow(non_normal_list_element) - 1) / sqrt(nrow(non_normal_list_element))


    ## Compare to critical value
    
    # An if statement to control on which tail the outlier was identified
    if (left_tail == TRUE) {
        message("Identifying the lowest outlier on the left tail...")
        
        # An if statement to determine if there is an outlier
        if (min_z_score > crit_z) {
            
            # Assign the outlier identity and index
            min_outlier_ident <- non_normal_list_element$f_constants[z_scores %in% min_z_score]
            min_outlier_index <- match(x = min_outlier_ident, non_normal_list_element$f_constants)

            message("The following minimum value was found to be an outlier: ", min_outlier_ident, "at the following row index: ", min_outlier_index, " and therefore will be removed.")
            
            # Assign the modified dataframe
            output_df <- non_normal_list_element[!non_normal_list_element$f_constants %in% min_outlier_ident, ]

            # Update the while loop condition
            while_loop_run <- TRUE

            
        } else {
            message("No lower outlier have been identified.")

            # Update the while loop condition
            while_loop_run <- FALSE
        }


    } else {
        message("Identifying the highest outlier on the right tail...")
        
        # An if statement to determine if there is an outlier
        if (max_z_score > crit_z) {
            
            # Assign the outlier identity and index
            max_outlier_ident <- non_normal_list_element$f_constants[z_scores %in% max_z_score]
            max_outlier_index <- match(x = max_outlier_ident, non_normal_list_element$f_constants)

            message("The following maximum value was found to be an outlier: ", max_outlier_ident, "\n", "at the following row index: ", max_outlier_index, " and therefore will be removed.")
            
            # Assign the modified dataframe
            output_df <- non_normal_list_element[!non_normal_list_element$f_constants %in% max_outlier_ident, ]

            # Update the while loop condition
            while_loop_run <- TRUE

            
        } else {
            message("No upper outlier have been identified.")

            # Update the while loop condition
            while_loop_run <- FALSE
        }


    }

    # Assign the output values to the output list
    output_list[[1]] <- output_df
    output_list[[2]] <- while_loop_run

    # Return the output list
    return(list(output_df, while_loop_run))

}




## A parametric outlier test function using the Grubbs test
## The function identifies adn removes outliers from the input df based on the f-constants
remove_outlier_f_const_param = function(normal_list_element, left_tail = TRUE) {
    
    
    ## Declare function variables

    # Define f-constants vector
    f_constants <- vector(mode = "numeric", length = nrow(normal_list_element))

    # Define an output dataframe
    output_df <- data.frame(matrix(nrow = nrow(normal_list_element), ncol = ncol(normal_list_element)))

    # Define the output list
    output_list <- vector(mode = "list", length = 2)

    # Grubb's test - lower outlier result
    grubbs_lower_res <- vector(mode = "list", length = 5)

    # Lower outlier identity
    lower_outlier_ident <- vector(mode = "numeric", length = 1)

    # Lower outlier index
    lower_outlier_index <- vector(mode = "numeric", length = 1)

    # Grubb's test - upper outlier result
    grubbs_upper_res <- vector(mode = "list", length = 5)

    # Upper outlier identity
    upper_outlier_ident <- vector(mode = "numeric", length = 1)

    # Upper outlier index
    upper_outlier_index <- vector(mode = "numeric", length = 1)

    # Variable to track if no outliers are found
    while_loop_run <- TRUE


    ## Calculate and assign the appropriate variables

    # Assign the f-constants to a new vector
    f_constants <- normal_list_element$f_constants

    #Assign the input dataframe as an output dataframe so in case of no outliers the original df will be returned
    output_df <- normal_list_element
    

    # This if statement controls if the upper or lower outlier is identified and removed
    if (left_tail == TRUE) {
        
        # Assign the value of Grubbs test on the left tail
        grubbs_lower_res <- outliers::grubbs.test(f_constants, opposite = TRUE)

        # This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_lower_res$p.value < 0.05) {

            # Assign outlier identity and index
            lower_outlier_ident <- as.numeric(stringr::str_extract(grubbs_lower_res$alternative, "[0-9]+\\.[0-9]+"))
            lower_outlier_index <- match(x = lower_outlier_ident, f_constants)

            # An if statement to make sure that the output is only returned if there was an outlier    
            if (length(lower_outlier_index) > 0) {
                message("The following lower outliers will be removed: ", lower_outlier_ident)
                
                # Assign the modified dataframe
                output_df <- normal_list_element[!f_constants %in% lower_outlier_ident, ]

                # Update the while loop condition
                while_loop_run <- TRUE
                
            }


        } else {
            message("No significant outlier found on the left tail.")
            
            # Update the while loop condition
            while_loop_run <- FALSE
            

        }


    } else {
        
        # Assign the value of Grubbs test on the right tail
        grubbs_upper_res <- outliers::grubbs.test(f_constants)

        # This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_upper_res$p.value < 0.05) {
            
            # Assign outlier identity and index
            upper_outlier_ident <- as.numeric(stringr::str_extract(grubbs_upper_res$alternative, "[0-9]+\\.[0-9]+"))
            upper_outlier_index <- match(x = upper_outlier_ident, f_constants)
                
            # An if statement to make sure that the output is only returned if there was an outlier
            if (length(upper_outlier_index) > 0) {
                message("The following upper outliers will be removed: ", upper_outlier_ident)
                
                # Assign the modified dataframe
                output_df <- normal_list_element[!f_constants %in% upper_outlier_ident, ]

                # Update the while loop condition
                while_loop_run <- TRUE
                
            }


        } else {
            message("No significant outlier found on the right tail.")
            
            # Update the while loop condition
            while_loop_run <- FALSE
            
        }
    }


    # Assign the output values to the output list
    output_list[[1]] <- output_df
    output_list[[2]] <- while_loop_run

    # Return the output list
    return(list(output_df, while_loop_run))
    
}




## This a wrapper function which binds the 3 outlier detection and removal functions together and runs them ina  while loop until there are no more
## outliers
## IMPORTANT NOTE: this function seems to work, but is a bit cobbled togeather as I'm not very familiar with while loops. Therefore please if you have more
## experience with them check and debug this bad boy! :)
outlier_cleaner = function(calculate_f_constants_output_list, is_data_normal_output_list, nonparam_test = "numeric_outlier_test") {
    for (index in seq_along(is_data_normal_output_list)) {
        if (is_data_normal_output_list[[index]]$p.value < 0.05) {
            message("It seems your f-constants show a non-normal distribution. The chosen non-parametric test will be used for outlier detection and removal.")

            if (nonparam_test == "numeric_outlier_test") {
                for (e in seq_along(calculate_f_constants_output_list)) {
                    calculate_f_constants_output_list[[e]] <- remove_outlier_f_const_nonparam(calculate_f_constants_output_list[[e]])

                }

                return(calculate_f_constants_output_list)
                
            } else {
                while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- modified_Zscore_test(calculate_f_constants_output_list[[i]], left_tail = TRUE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

                while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- modified_Zscore_test(calculate_f_constants_output_list[[i]], left_tail = FALSE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

                return(calculate_f_constants_output_list)
            }


        } else {
            message("It seems your f-constants show a normal distribution. A parametric Grubbs test will be used for outlier detection and removal.")

            while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- remove_outlier_f_const_param(calculate_f_constants_output_list[[i]], left_tail = TRUE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

                while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- remove_outlier_f_const_param(calculate_f_constants_output_list[[i]], left_tail = FALSE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

            return(calculate_f_constants_output_list)
        }
    }
}




for (e in seq_along(f2)) {
                    f3[[e]] <- remove_outlier_f_const_nonparam(f3[[e]])

}











remove_outlier_f_const_param = function(normal_list_element, left_tail = TRUE) {
    
    
    ## Declare function variables

    # Define f-constants vector
    f_constants <- vector(mode = "numeric", length = nrow(normal_list_element))

    # Define an output dataframe
    output_df <- data.frame(matrix(nrow = nrow(normal_list_element), ncol = ncol(normal_list_element)))

    # Grubb's test - lower outlier result
    grubbs_lower_res <- vector(mode = "list", length = 5)

    # Lower outlier identity
    lower_outlier_ident <- vector(mode = "numeric", length = 1)

    # Lower outlier index
    lower_outlier_index <- vector(mode = "numeric", length = 1)

    # Grubb's test - upper outlier result
    grubbs_upper_res <- vector(mode = "list", length = 5)

    # Upper outlier identity
    upper_outlier_ident <- vector(mode = "numeric", length = 1)

    # Upper outlier index
    upper_outlier_index <- vector(mode = "numeric", length = 1)

    # Variable to track if no outliers are found
    while_loop_run <- TRUE


    ## Calculate and assign the appropriate variables

    # Assign the f-constants to a new vector
    f_constants <- normal_list_element$f_constants

    #Assign the input dataframe as an output dataframe so in case of no outliers the original df will be returned
    output_df <- normal_list_element
    

    # This if statement controls if the upper or lower outlier is identified and removed
    if (left_tail == TRUE) {
        
        # Assign the value of Grubbs test on the left tail
        grubbs_lower_res <- outliers::grubbs.test(f_constants, opposite = TRUE)

        # This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_lower_res$p.value < 0.05) {

            # Assign outlier identity and index
            lower_outlier_ident <- as.numeric(stringr::str_extract(grubbs_lower_res$alternative, "[0-9]+\\.[0-9]+"))
            lower_outlier_index <- match(x = lower_outlier_ident, f_constants)

            # An if statement to make sure that the output is only returned if there was an outlier    
            if (length(lower_outlier_index) > 0) {
                message("The following lower outliers will be removed: ", lower_outlier_ident)
                
                # Assign the modified dataframe
                output_df <- normal_list_element[!f_constants %in% lower_outlier_ident, ]
                print(output_df)

                # Update the while loop condition
                while_loop_run <- TRUE
                
            }


        } else {
            message("No significant outlier found on the left tail.")
            
            # Update the while loop condition
            while_loop_run <- FALSE
            

        }


    } else {
        
        # Assign the value of Grubbs test on the right tail
        grubbs_upper_res <- outliers::grubbs.test(f_constants)

        # This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_upper_res$p.value < 0.05) {
            
            # Assign outlier identity and index
            upper_outlier_ident <- as.numeric(stringr::str_extract(grubbs_upper_res$alternative, "[0-9]+\\.[0-9]+"))
            upper_outlier_index <- match(x = upper_outlier_ident, f_constants)
                
            # An if statement to make sure that the output is only returned if there was an outlier
            if (length(upper_outlier_index) > 0) {
                message("The following upper outliers will be removed: ", upper_outlier_ident)
                
                # Assign the modified dataframe
                output_df <- normal_list_element[!f_constants %in% upper_outlier_ident, ]
                print(output_df)

                # Update the while loop condition
                while_loop_run <- TRUE
                
            }


        } else {
            message("No significant outlier found on the right tail.")
            
            # Update the while loop condition
            while_loop_run <- FALSE
            
        }
    }

    return(list(output_df, while_loop_run))
    
}



while_loop_run <- TRUE
while (while_loop_run) {
    for (i in seq_along(f3)) {
        result <- remove_outlier_f_const_param(f3[[i]], left_tail = FALSE)
        f3[[i]] <- result[[1]]
        print(result)
        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
    }
}


while_loop_run <- TRUE
while (while_loop_run) {
    for (i in seq_along(f3)) {
        result <- modified_Zscore_test(f3[[i]], left_tail = FALSE)
        f3[[i]] <- result[[1]]
        print(result)
        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
    }
}




## Remove the outliers using Grubb's test (part of the outlier package) (parametric test!)
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





## Non-parametric outlier test - Numeric outlier test
remove_outlier_f_const_nonparam = function(calculate_f_constants_output_list) {
    ## Declare a list to modify
    list_to_modify <- calculate_f_constants_output_list

    # Declare the output list
    f_constant_outlier_free_measurements <- vector(mode = "list", length = length(list_to_modify))

    ## Carry out the calculations, outlier identification and removal for each input list element
    for (list_i in seq_along(list_to_modify)) {
        ## Declare dynamic variables

        # Declare the quartile variables
        Q1 <- vector(mode = "numeric", length = 1)
        Q3 <- vector(mode = "numeric", length = 1)

        # Declare a vector for the IQR - inter quartile range
        IQR <- vector(mode = "numeric", length = 1)

        # Declare vectors containing the upper and lower boundaries
        upper_boundary <- vector(mode = "numeric", length = 1)
        lower_boundary <- vector(mode = "numeric", length = 1)

        # Declare vectors for the upper and lower outliers
        upper_outlier <- vector(mode = "numeric")
        lower_outlier <- vector(mode = "numeric")

        # Order the constants in am ascending order
        ordered_constants <- list_to_modify[[list_i]]$f_constants[order(list_to_modify[[list_i]]$f_constants)]


        # Calculate the quartile ranges and allocate them to the appropriate list elements
        Q1 <- quantile(ordered_constants, 0.25, na.rm = TRUE)
        Q3 <- quantile(ordered_constants, 0.75, na.rm = TRUE)

        # Calculate the IQR 
        IQR <- Q3 - Q1

        # Calculate the upper and lower boundaries
        upper_boundary <- Q3 + (1.5 * IQR)
        lower_boundary <- Q1 - (1.5 * IQR)

        if (lower_boundary < 0) {
            lower_boundary <- 0
        }

        # Identify the upper outliers
        upper_outlier <- ordered_constants[ordered_constants > upper_boundary]

        # Identify the lower outliers
        lower_outlier <- ordered_constants[ordered_constants < lower_boundary]
        
        # An if statement to check if there were any upper outliers and if yes to remove them
        if (length(upper_outlier) == 0) {
            message("No outlier found on the right (upper) tail for input list element: ", list_i)
        } else {
            message("The following elements are outliers on the right (upper) tail, and will be removed: \n", upper_outlier)
            
            #Remove the upper outliers
            list_to_modify[[list_i]] <- list_to_modify[[list_i]][!list_to_modify[[list_i]]$f_constants %in% upper_outlier, ]
        }

        # An if statement to check if there were any lower outliers and if yes to remove them
        if (length(lower_outlier) == 0) {
            message("No outlier found on the left (lower) tail for input list element: ", list_i)
        } else {
            message("The following elements are outliers on the left (lower) tail, and will be removed: \n", lower_outlier)

            #Remove the lower outliers
            list_to_modify[[list_i]] <- list_to_modify[[list_i]][!list_to_modify[[list_i]]$f_constants %in% lower_outlier, ]
        }

    }

    # Assign the f-constant outlier free list element to the output list
    f_constant_outlier_free_measurements <- list_to_modify

    # Return the output list
    return(f_constant_outlier_free_measurements)
    
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
