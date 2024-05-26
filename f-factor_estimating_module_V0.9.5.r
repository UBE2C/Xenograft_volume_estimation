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

            } 
        }

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

# Define the preferred input directory path
pref_input_path <- paste0(script_dir_path, "/", "Input_files")

# Define the preferred path to the intermediate input/output file directory 
pref_output_path <- paste0(script_dir_path, "/", "Intermediate_IO")

# Define the preferred final output directory path  for later use
final_output_path <- paste0(script_dir_path, "/", "Final_output_files")

# Directory management function
# This function will look for the preferred input and output libraries and will create them if they are missing upon user request
dir_management = function(input_1 = NULL, input_2 = NULL) {
    ##Define the function variables
    
    #Declare a vector containing the preferred directory paths
    dir_paths <- c(pref_input_path, pref_output_path, final_output_path)
    
    #This for loop will traverse the directory paths vector to feed the paths to the if statement below
    for (path in dir_paths) {
        
        #The main if statement which will check if a directory exists
        if (!dir.exists(path)) {
            message("The following directory does not exists but is recommended for better organization: ", path, "\n",
                "would you like to create it? \n",
                "1: Yes \n",
                "2: No \n",
                "(select a number)")
            
            #The scan function will stop execution and will wait for user input, then it will set the value of
            #the input_1 variable based on user input
            input_1 <- scan(file = "stdin", what = numeric(), n = 1, quiet = TRUE)

            #An inner if statement to control the next steps based on user input
            if (input_1 == 1) {
                dir.create(path)

                message("The requested directories were successfully created!")
            } else {
                message("The directory ", path, "was not requested. \n",
                    "The default paths will be used")
            }
        }
    }

    #An if statement to offer the user to re-organize the input files
    if (!is.null(input_1) && input_1 == 1) {
        message("Now that you created the new directories would like to stop and place the input files in the proper directories? \n",
            "1: Yes \n",
            "2: No \n",
            "(select a number)")

        input_2 <- scan(file = "stdin", what = numeric(), n = 1, quiet = TRUE)
    }

    #An if statement to stop execution if the user choses so
    if (!is.null(input_2) && input_2 == 1) {
        quit(save = "no")

    }

    #An if statement to set the the default input and output paths to the preferred ones if they exists
    if (dir.exists(dir_paths[1]) == TRUE && dir.exists(dir_paths[2]) == TRUE) {
        # Define the preferred input directory path
        assign("def_input_path", pref_input_path, envir = parent.frame())

        # Define the final output directory path
        assign("def_output_path", pref_output_path, envir = parent.frame())

    } else {
        # Define the input directory path
        assign("def_input_path", script_dir_path, envir = parent.frame())

        # Define the final output directory path
        assign("def_output_path", script_dir_path, envir = parent.frame())

    }
    
}


################################################# Section end #################################################






################################################# MARK: call the Dir_management function #################################################
                                                #                                        #




## Calling the directory management function 
## NOTE: tis function must be called before the option parsing, so it can properly assign and modify the input and output directories
## before the options parsing happens!
dir_management()


################################################# Section end #################################################






################################################# MARK: CLI arguments #################################################
                                                #                     #


## Define the default output file names
def_output_name <- "not_defined_yet"



## Add command line arguments
options_list <- list(
    optparse::make_option(opt_str = c("-p", "--input_path"), action = "store", type = "character", default = def_input_path,
    help = "This argument takes a character string which defines the path to the input .csv files ready to be processed. \n
    By default the script sets the local directory as path."),

    optparse::make_option(opt_str = c("-o", "--output_path"), action = "store", type = "character", default = def_output_path,
    help = "This argument takes a character string which defines the output path to the output .csv files after processing. \n
    By default the path is called Intermediate_IO, which is a directory found under the local directory."),

    optparse::make_option(opt_str = c("-n", "--output_name"), action = "store", type = "character", default = def_output_name,
    help = "This argument takes a character string which defines the name of the output .csv files. \n
    By default the script sets the name to processed_output_table."),

    optparse::make_option(opt_str = c("-v", "--verbose"), action = "store_true", default = FALSE, dest = "verbose",
    help = "This argument controls if the program returns additional information like the outputs of the sub-functions, messages and warnings. \n
    By default the option is set to FALSE."),

    optparse::make_option(opt_str = c("-q", "--quiet"), action = "store_true", default = FALSE, dest = "quiet"),

    optparse::make_option(opt_str = c("-np", "--nonparam_test"), default = "numeric_outlier_test"),

    optparse::make_option(opt_str = c("-ol", "--outlier_handling"), default = "detect"),

    optparse::make_option(opt_str = c("-cr", "--correction"))

)

# Create a program description
prog_descr <- c(paste0(
                "For a standalone run use the 'Rscript COMIX-data_processer_V[option] [optional args]' syntax."))

# Parse the arguments to the arguments variable
arguments <- optparse::parse_args(object = optparse::OptionParser(option_list = options_list, 
                                    description = prog_descr),
                                    args = commandArgs(trailingOnly = TRUE),
                                    print_help_and_exit = TRUE,
                                    positional_arguments = FALSE,
                                    convert_hyphens_to_underscores = FALSE)


################################################# Section end #################################################





################################################# MARK: Data reading #################################################
                                                #   and management   #




## Loading and cleaning the uCT data files needed for the calculations
## NOTE: the first sub-function needs a new file with the uCT measurements 


# This function loads the uCT data from the intermediate I/O folder for further processing
read_uCT_data = function(data_path, separator) {
    ##Declare function variables
    
    #Declare a list onto which the uCT data will be added
    uCT_measurements <- list()
    

    ##Start reading  the appropriate files

    #List the files found in the intermediate I/O folder
    intermediate_IO_files <- list.files(path = data_path)
    
    #Grep the uCT file names for the checking if statement
    uCT_file_names <- intermediate_IO_files[grep(pattern = "uCT", x = intermediate_IO_files, ignore.case = TRUE)]

    #This main if statement will check if the correct files are in the folder
    if (length(uCT_file_names) > 0) {
        
        #Read each uCT .csv file and add them to the output list
        for (item in intermediate_IO_files) {
            if (grepl(pattern = "uCT", x = item, ignore.case = TRUE) == TRUE) {
                uCT_measurements[[item]] <- read.csv(file = paste0(data_path, "/", item), sep = separator)
            }
        }

        #Name the return list elements based on the original file names
        names(uCT_measurements) <- uCT_file_names

        
        ##Return the uCT data
        return(uCT_measurements)

    } else {
        warning("No uCT input table was found in the ", data_path, " directory.")
    }

}


# This function cleans the uCT data by removing entries without corresponding caliper measurements
clean_uCT_data = function(read_uCT_output_list, verbose) {
    ##Declare function variables
    
    #Declare a new list for the clean uCT dataframes
    clean_uCT_measurement <- list()

    #Declare an intermediate list to modify
    temp_uCT_df <- data.frame()

    #Define a vector containing the bools for the columns which contains only NAs if any
    only_na_cols <- vector(mode = "logical")
    

    ##Clean the dataframes by removing rows containing sparse NA values

    #This main for loop will walk through the uCT list of dataframes and removes rows where there are sparse NA values       
    for (index in seq_along(read_uCT_output_list)) {
        #Initialize an intermediate temporary dataframe to store the list elements one-by-one
        temp_uCT_df <- read_uCT_output_list[[index]]

        #Loop through the measurement columns and look for columns where are some measurements with NAs mixed in
        for (col in seq(from = 4, to = ncol(temp_uCT_df), by = 1)) {
            #This if statement will look for columns where are some measurements with NAs mixed in (sparse NAs)
            if (length(unique(is.na(read_uCT_output_list[[index]][, col]))) > 1) {
                #A message to warn the user that NA's were found and the appropriate sample row will be removed
                if (verbose == TRUE) {
                    warning("Some missing values were found in the following table: \n", names(read_uCT_output_list)[index], "\n", "at the following column: \n", colnames(temp_uCT_df)[col], "\n",
                    "at the following sample: \n", temp_uCT_df$Mouse_ID[is.na(temp_uCT_df[, col])], "\n",
                    "The affected sample will be removed for the analysis.")

                }
                
                #Remove the appropriate sample row and update the temp dataframe with the new state
                temp_uCT_df <- temp_uCT_df[!is.na(temp_uCT_df[, col]), ]  

            } 

        }

        #This for loop will loop through the columns of the temp_datafame and flags each columns which contains only NAs for removal
        for (i in seq_len(ncol(temp_uCT_df))) {
            only_na_cols[i] <- length(unique(is.na(temp_uCT_df[, i]))) == 1 && unique(is.na(temp_uCT_df[, i])) == TRUE

        }

        #Remove the columns with only NA values
        if (verbose == TRUE) {
            warning("In the following table: \n", names(read_uCT_output_list)[index], "\n", "the following columns: \n", colnames(temp_uCT_df)[only_na_cols], "\n",
                    "contain only missing values therefore they will be removed for the analysis")
        }
        
        temp_uCT_df <- temp_uCT_df[, !only_na_cols]
        

        #Reset the row numbers for the modified dataframe
        rownames(temp_uCT_df) <- NULL

        #Assign the modified dataframe to the clean list
        clean_uCT_measurement[[index]] <- temp_uCT_df

    }

    #Name the elements of the return list based on the elements of the input list
    names(clean_uCT_measurement) <- names(read_uCT_output_list)


    ##Return the clean uCT output list
    return(clean_uCT_measurement)

}




## Loading the processed caliper measurement data files needed for the calculations
## NOTE: this sub-function needs the clean output from the previous module


# This function loads the cleaned output files from the Data-processer module
read_caliper_data = function(data_path, separator) {
    ##Declare function variables
    
    #Declare a list onto which the caliper measurements will be added
    caliper_measurements <- list()


    ##Start reading the appropriate files

    #List the files found in the intermediate I/O folder
    intermediate_IO_files <- list.files(path = data_path)
    
    #Grep the uCT file names for the checking if statement
    processed_caliper_file_names <- intermediate_IO_files[grep(pattern = "caliper", x = intermediate_IO_files, ignore.case = TRUE)]
    
    #This main if statement will check if the correct files are in the folder
    if (length(processed_caliper_file_names) > 0) {
        #Add the measurements to the list with a for loop
        for (item in intermediate_IO_files) {
            if (grepl("caliper", item, ignore.case = TRUE) == TRUE) {
                caliper_measurements[[item]] <- read.csv(file = paste0(data_path, "/", item), sep = separator)

            }
            
        }
    } else {
        warning("No processed caliper measurement files was found in the ", data_path, " directory.")

    }

    #Name the elements of the caliper_measurements list based on the original processed I/O files
    names(caliper_measurements) <- processed_caliper_file_names


    ##Return the caliper_measurement list
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

    #Initialize a list onto which the row trimmed down caliper measurements will be added
    trimmed_caliper_list <- vector(mode = "list")
    
    #Initialize temporary dataframes on which the processing can take place
    temp_uct_df <- data.frame()
    temp_calip_df <- data.frame()

    #Initialize a vector which will hold the name of the caliper list element which matches the uCT list element
    list_element_name <- vector(mode = "character", length = 1)    


    ##Trim the caliper measurements

    #The main body of the function responsible for the trimming of the caliper list elements according to the uCT list elements
    #The outer for loop traverses the uCT measurement list 
    for (element in seq_along(clean_uCT_data_output_list)) {
        #Assign the list elements (dataframes) into the appropriate temporary dataframe variable
        temp_uct_df <- clean_uCT_data_output_list[[element]]

        #Initialize a list onto which the sampling dates will be added (NOTE: the function assumes that the sampling dates between the uCT and caliper measurements are the same!)
        uCT_sampling_dates <- unlist(stringr::str_extract_all(string = colnames(temp_uct_df), pattern = "_[0-9\\.]+"))
        uCT_sampling_dates <<- unlist(stringr::str_extract_all(string = colnames(temp_uct_df), pattern = "_[0-9\\.]+"))

        
        #This inner loop traverses the caliper measurements list
        for (item in seq_along(read_caliper_data_output_list)) {
            #This if statement is responsible for matching the current element of the uCT list with the corresponding element of the caliper list
            if (unique(temp_uct_df$Treatment_group_id) == unique(read_caliper_data_output_list[[item]]$Treatment_group_id)) {
                #Assign the corresponding caliper list element (dataframes) into the appropriate temporary dataframe variable
                temp_calip_df <- read_caliper_data_output_list[[item]]

                #Assign the corresponding caliper list element name into the temporary list element name variable
                list_element_name <- names(read_caliper_data_output_list)[item]
            }

            #Trim the corresponding caliper measurement dataframe by the correct rows
            temp_calip_df <- temp_calip_df[temp_calip_df$Mouse_ID %in% temp_uct_df$Mouse_ID, ]

            #Initialize a dynamic variable to store the positions of the corresponding uCT dates (using the trimmed uCT df) 
            matching_cols <- vector(mode = "numeric", length = length(uCT_sampling_dates))

            #This for loop will pull the positions of the matching columns based on the dates of the uCT and caliper measurements
            for (i in seq_along(uCT_sampling_dates)){
                #Assign the positions of the matching columns
                matching_cols[i] <- grep(pattern = uCT_sampling_dates[i], x = colnames(temp_calip_df))
            }
            
            #Trim the corresponding caliper measurement dataframe by the correct columns
            temp_calip_df <- temp_calip_df[, c(1, 2, 3, matching_cols)]

        }

        #Assign the trimmed caliper measurement to the output list
        trimmed_caliper_list[[element]] <- temp_calip_df

        #Add the matching name to each list element
        names(trimmed_caliper_list)[element] <- list_element_name

    }


    ##Return the the fully trimmed data as a list
    return(trimmed_caliper_list)

}


#Bind the trimmed caliper measurement dataframes together by columns and create a clean unified df
#This function will carry out the column bind and organizes the dataframes according to the uCT measurements list, then creates a list of unified dataframes
#using the selected column of the cuCT measurements and the transposed bound caliper measurements
bind_and_unify_measurements = function(fit_caliper_measurements_output_list, clean_uCT_data_output_list) {
    
    
    ##Define the variables used in the function

    #Initialize a unified dataframe list, which will unify the bound caliper measurements and selected columns from the clean_uCT_list
    unified_df_list <- list()


    ##Start the processing steps

    
    #Unify the bound caliper measurements and selected columns from the clean_uCT_list
    for (i in seq_along(fit_caliper_measurements_output_list)) {
        unified_df_list[[i]] <- cbind(fit_caliper_measurements_output_list[[i]], clean_uCT_data_output_list[[i]][, c(grep(pattern = "volume", x = colnames(clean_uCT_data_output_list[[i]]), ignore.case = TRUE))])
    }

    #Name the output list elements
    names(unified_df_list) <- names(fit_caliper_measurements_output_list)

    ##Return the bound_df_list
    return(unified_df_list)

}


################################################# Section end #################################################





################################################    MARK: f-constant    #################################################
                                                #  calculation and eval #




## Calculate the f-constants to each uCT-measurement-caliper measurement group


# NOTE: original formula V = (pi/6) * f * (l * w)^(3/2)
# NOTE: formula solved to f f = V / (pi/6) * (l * w)^(3/2)


# This function calculates the sample mean f-constants based on the tumor volumes form the uCTs and the caliper measurements from the LxW values
calculate_f_constants = function(bind_and_unify_measurements_output_list) {
    ##Define the variables used in the function

    #Initialize a list which will hold the calculated f-constants for each sample of each dataframe
    f_constants <- vector(mode = "list", length = length(bind_and_unify_measurements_output_list))

    #Initialize a new list which will hold the new column names for the calculated f-constants (designated by date)
    new_col_names <- vector(mode = "list", length = length(bind_and_unify_measurements_output_list))

    #Initialize the final output list
    output_list <- vector(mode = "list", length = length(bind_and_unify_measurements_output_list))
    

    ##Separate the measurements by date and calculate the corresponding f-constant

    #This loop will traverse the input list of dataframes and splits down the input list into individual dataframes
    for (element in seq_along(bind_and_unify_measurements_output_list)) {
        ##Define dynamic variables

        #Initialize a dynamic temp_dataframe to store currently processed list element
        temp_unif_df <- bind_and_unify_measurements_output_list[[element]]

        #Initialize a dynamic temp_dataframe to store currently processed list element
        temp_trim_df <- data.frame()

        #Initialize vectors for the the single LxW and volume
        temp_LxW_values <- vector(mode = "numeric", length = nrow(temp_unif_df))
        temp_ct_volumes <- vector(mode = "numeric", length = nrow(temp_unif_df))

        #Initialize a dynamic date variable to store the unique dates found in the current df
        sampling_dates <- unique(unlist(stringr::str_extract_all(string = colnames(temp_unif_df), pattern = "_[0-9\\.]+")))

        #The inner loop will traverse the extracted sampling dates to split the dataframes down to pairs of LxV-volume columns for processing
        for (item in seq_along(sampling_dates)) {
            #Split down the current measurement dataframe into a single pair of LxW - Volume columns which belongs to the same sampling date
            temp_trim_df <- temp_unif_df[, grep(pattern = sampling_dates[item], x = colnames(temp_unif_df))]

            #Further split the previously created trimmed dataframes into a vector of LxW values
            temp_LxW_values <- temp_trim_df[, grep(pattern = "LxW", x = colnames(temp_trim_df), ignore.case = TRUE)]

            #Further split the previously created trimmed dataframes into a vector of Volume values
            temp_ct_volumes <- temp_trim_df[, grep(pattern = "volume", x = colnames(temp_trim_df), ignore.case = TRUE)]

            #Use the LxW and Volume vectors to calculate the f-constant for each sample and assign the results to the f-constants list (creates a nested list)
            f_constants[[element]][[item]] <- temp_ct_volumes / ((pi / 6) * (temp_LxW_values)^(3 / 2))

            
        }

        #Assign the new f-constant column names (based on sampling dates) to the initialized variable
        new_col_names[[element]] <- paste0("f_const", sampling_dates)
        
    }

    #This loop traverses the f-constants nested list to take each individual f-factor vector and add it to the input list dataframes as columns
    for (element in seq_along(f_constants)) {
        ##Define dynamic variables

        #Initialize a dynamic dataframe holding the current input list element (dataframe)
        temp_modif_df <- bind_and_unify_measurements_output_list[[element]]
        
        #This inner loop traverses the nested items in the f-constants list (the calculated f-constants based on sampling dates)
        for (item in seq_along(f_constants[[element]])) {
            #Bind the f-constant vector to the dynamic dataframe 
            temp_modif_df <- cbind(temp_modif_df, f_constants[[element]][[item]])

        }
        
        #Assign the appended dataframes as the appropriate output list element
        output_list[[element]] <- temp_modif_df

        #Change the names of the new columns in each list element(dataframe) to the appropriate name
        colnames(output_list[[element]])[grep(pattern = "f_constants", x = colnames(output_list[[element]]))] <- new_col_names[[element]]
    }

    #This loop traverses the output list to calculate a mean of the f-constants calculated to each date
    for (element in seq_along(output_list)) {
        ##Define dynamic variables

        #Initialize a dynamic dataframe which will hold the mean f-constants
        temp_element_df <- output_list[[element]]

        #This inner loop traverses the rows of every dataframe in order to trim them down to the calculated f-constants columns and calculate the mean
        #f-constant values for each sample
        for (row in seq_len(nrow(temp_element_df))) {
            #Split the dataframes to the f-constant columns and calculate the mean values for each sample (row)
            temp_element_df$f_const_mean[row] <- mean(as.numeric(temp_element_df[row, grep(pattern = "f_const", x = colnames(temp_element_df))]))

        }

        #Re-assign the modified dataframes holding the mean_f-constants to the final output list
        output_list[[element]] <- temp_element_df
        
    }

    #Assign the original list element names to the output list
    names(output_list) <- names(bind_and_unify_measurements_output_list)


    ##Return the final output list
    return(output_list)

}


################################################# Section end #################################################





################################################      MARK: Outlier calculation      #################################################
                                                #  and removal among the f-constants #




## Determine if the data is normally distributed
is_data_normal = function(calculate_f_constants_output_list) {
    

    ##Declare dynamic variables
    shapiro_results <- vector(mode = "list", length = length(calculate_f_constants_output_list))


    ##Do the normality test

    #This for loop will carry out a Shapiro-Wilk normality test
    for (i in seq_along(calculate_f_constants_output_list)) {
        shapiro_results[[i]] <- shapiro.test(calculate_f_constants_output_list[[i]]$f_const_mean)
    }


    ##Return the shapiro test results list
    return(shapiro_results)

}




## Non-parametric outlier test - Numeric outlier test (this seems to be the more safe alternative as it is a standard test)
## The function identifies and removes outliers from the input df based on the f-constants
remove_outlier_f_const_NOTest = function(non_normal_list_element) {
    ##Declare a list to modify
    list_element_to_modify <- non_normal_list_element

    #Declare the output list
    f_constant_outlier_free_measurements <- data.frame(matrix(nrow = nrow(list_element_to_modify), ncol = ncol(list_element_to_modify)))

    ##Declare variables

    #Declare the quartile variables
    Q1 <- vector(mode = "numeric", length = 1)
    Q3 <- vector(mode = "numeric", length = 1)

    #Declare a vector for the IQR - inter quartile range
    IQR <- vector(mode = "numeric", length = 1)

    #Declare vectors containing the upper and lower boundaries
    upper_boundary <- vector(mode = "numeric", length = 1)
    lower_boundary <- vector(mode = "numeric", length = 1)

    #Declare vectors for the upper and lower outliers
    upper_outlier <- vector(mode = "numeric")
    lower_outlier <- vector(mode = "numeric")

    #Order the constants in am ascending order
    ordered_constants <- list_element_to_modify$f_const_mean[order(list_element_to_modify$f_const_mean)]


    #Calculate the quartile ranges and allocate them to the appropriate list elements
    Q1 <- quantile(ordered_constants, 0.25, na.rm = TRUE)
    Q3 <- quantile(ordered_constants, 0.75, na.rm = TRUE)

    #Calculate the IQR 
    IQR <- Q3 - Q1

    #Calculate the upper and lower boundaries
    upper_boundary <- Q3 + (1.5 * IQR)
    lower_boundary <- Q1 - (1.5 * IQR)

    if (lower_boundary < 0) {
        lower_boundary <- 0
    }

    #Identify the upper outliers
    upper_outlier <- ordered_constants[ordered_constants > upper_boundary]

    #Identify the lower outliers
    lower_outlier <- ordered_constants[ordered_constants < lower_boundary]
        
    #An if statement to check if there were any upper outliers and if yes to remove them
    if (length(upper_outlier) == 0) {
        message("No outlier found on the right (upper) tail for input list element. \n")
    } else {
        message("The following elements are outliers on the right (upper) tail, and will be removed: \n", upper_outlier, "\n")
            
        #Remove the upper outliers
        list_element_to_modify <- list_element_to_modify[!list_element_to_modify$f_const_mean %in% upper_outlier, ]
    }

    #An if statement to check if there were any lower outliers and if yes to remove them
    if (length(lower_outlier) == 0) {
        message("No outlier found on the left (lower) tail for input list element. \n")
    } else {
        message("The following elements are outliers on the left (lower) tail, and will be removed: \n", lower_outlier, "\n")

        #Remove the lower outliers
        list_element_to_modify <- list_element_to_modify[!list_element_to_modify$f_const_mean %in% lower_outlier, ]
    }



    #Assign the f-constant outlier free list element to the output list
    f_constant_outlier_free_measurements <- list_element_to_modify


    ##Return the output list
    return(f_constant_outlier_free_measurements)
    
}




## Non-parametric outlier test - Numeric outlier test (NOT) (this seems to be the more safe alternative as it is a standard test)
## The function identifies outliers from the input df based on the f-constants
detect_outlier_f_const_NOTest = function(non_normal_list_element) {
    ##Declare the output list
    outliers_list <- list("upper_outliers" = NULL, "lower_outliers" = NULL)


    ##Declare variables

    #Declare the quartile variables
    Q1 <- vector(mode = "numeric", length = 1)
    Q3 <- vector(mode = "numeric", length = 1)

    #Declare a vector for the IQR - inter quartile range
    IQR <- vector(mode = "numeric", length = 1)

    #Declare vectors containing the upper and lower boundaries
    upper_boundary <- vector(mode = "numeric", length = 1)
    lower_boundary <- vector(mode = "numeric", length = 1)

    #Declare vectors for the upper and lower outliers
    upper_outlier <- vector(mode = "numeric")
    lower_outlier <- vector(mode = "numeric")

    #Order the constants in am ascending order
    ordered_constants <- non_normal_list_element$f_const_mean[order(non_normal_list_element$f_const_mean)]


    #Calculate the quartile ranges and allocate them to the appropriate list elements
    Q1 <- quantile(ordered_constants, 0.25, na.rm = TRUE)
    Q3 <- quantile(ordered_constants, 0.75, na.rm = TRUE)

    #Calculate the IQR 
    IQR <- Q3 - Q1

    #Calculate the upper and lower boundaries
    upper_boundary <- Q3 + (1.5 * IQR)
    lower_boundary <- Q1 - (1.5 * IQR)

    if (lower_boundary < 0) {
        lower_boundary <- 0
    }

    #Identify the upper outliers
    upper_outlier <- ordered_constants[ordered_constants > upper_boundary]

    #Identify the lower outliers
    lower_outlier <- ordered_constants[ordered_constants < lower_boundary]
        
    #An if statement to check if there were any upper outliers and if yes to remove them
    if (length(upper_outlier) == 0) {
        message("No outlier found on the right (upper) tail for input list element. \n")

    } else {
        message("The following elements are outliers on the right (upper) tail: \n", upper_outlier, "\n")

        #Assign the outlier f-constants to the output list
        outliers_list[[1]] <- upper_outlier
            
    }

    #An if statement to check if there were any lower outliers and if yes to remove them
    if (length(lower_outlier) == 0) {
        message("No outlier found on the left (lower) tail for input list element. \n")

    } else {
        message("The following elements are outliers on the left (lower) tail: \n", lower_outlier, "\n")

        #Assign the outlier f-constants to the output list
        outliers_list[[2]] <- lower_outlier
    }

    if (length(upper_outlier) != 0 || length(lower_outlier != 0)) {
        ##Return the output list
        return(outliers_list)

    }
    
    
}




## Modified z-score test (non-parametric) - this test is basically the Grubb's test, just with the median and MAD instead of mean and SD
## The ida was taken from the following publication - DOI: (https://doi.org/10.1515/dema-2021-0041), however the publication used the SD with the median
## and not the MAD. The threshold to which the z-score is compared is calculated with the simple N-1/sqrt(N) formula taken from Graphpad Prism's test description.
## link: (https://www.graphpad.com/guides/prism/latest/statistics/stat_detecting_outliers_with_grubbs.htm)
## The function identifies and removes outliers from the input df based on the f-constants
remove_outlier_f_const_mZscore_test = function(non_normal_list_element, left_tail = TRUE) {
    
    
    ##Define function variables

    #Define an output dataframe
    output_df <- data.frame(matrix(nrow = nrow(non_normal_list_element), ncol = ncol(non_normal_list_element)))

    #Define the output list
    output_list <- vector(mode = "list", length = 2)
    
    #Define the f_constant vector for sorting
    ordered_f_constants <- vector(mode = "numeric", length = nrow(non_normal_list_element))
    
    #Define the median f-constant
    f_const_median <- vector(mode = "numeric", length = 1)

    #Define the median absolute deviation (MAD)
    f_const_mad <- vector(mode = "numeric", length = 1)

    #Define a vector for the z-scores
    z_scores <- vector(mode = "numeric", length = nrow(non_normal_list_element))

    #Define a vector for the minimum z-score
    min_z_score <- vector(mode = "numeric", length = 1)

    #Define a vector for the maximum z-score
    max_z_score <- vector(mode = "numeric", length = 1)

    #Define the critical z-value
    crit_t <- vector(mode = "numeric", length = 1)

    #Max outlier identity
    upper_outlier_ident <- vector(mode = "numeric", length = 1)

    #Max outlier index
    max_outlier_index <- vector(mode = "numeric", length = 1)

    #Min outlier identity
    lower_outlier_ident <- vector(mode = "numeric", length = 1)

    #Min outlier index
    min_outlier_index <- vector(mode = "numeric", length = 1)

    #Define the degrees of freedom
    deg_of_f <- vector(mode = "numeric", length = 1)

    #Variable to track if no outliers are found
    while_loop_run <- TRUE


    ##Calculate the defined variables

    #Assign and sort the f_constants
    ordered_f_constants <- non_normal_list_element$f_const_mean[order(non_normal_list_element$f_const_mean)]

    #Assign the input dataframe as an output dataframe so in case of no outliers the original df will be returned
    output_df <- non_normal_list_element

    #Calculate the median of the f-constants
    f_const_median <- median(ordered_f_constants, na.rm = TRUE)

    #Calculate the MAD
    f_const_mad <- mad(x = non_normal_list_element$f_const_mean, na.rm = TRUE)

    #Calculate the minimum z-score
    min_z_score <- (f_const_median - min(non_normal_list_element$f_const_mean)) / f_const_mad

    #Calculate all the z-scores
    for (i in seq_along(non_normal_list_element$f_const_mean)) {
        z_scores[i] <- abs((non_normal_list_element$f_const_mean[i] - f_const_median) / f_const_mad)
    
    }

    #Calculate the maximum z-score
    max_z_score <- (max(non_normal_list_element$f_const_mean) - f_const_median) / f_const_mad

    #Calculate the degrees of freedom
    deg_of_f <- nrow(non_normal_list_element) - 1

    #Calculate the critical value of t
    crit_t <- qt(1 - 0.05 / 2, deg_of_f)


    ##Compare to critical value
    
    #An if statement to control on which tail the outlier was identified
    if (left_tail == TRUE) {
        message("Identifying the lowest outlier on the left tail...")
        
        #An if statement to determine if there is an outlier
        if (min_z_score > crit_t) {
            
            #Assign the outlier identity and index
            lower_outlier_ident <- non_normal_list_element$f_const_mean[z_scores %in% min_z_score]
            min_outlier_index <- match(x = lower_outlier_ident, non_normal_list_element$f_const_mean)

            message("The following minimum value was found to be an outlier: ", lower_outlier_ident, "at the following row index: ", min_outlier_index, " and therefore will be removed. \n")
            
            #Assign the modified dataframe
            output_df <- non_normal_list_element[!non_normal_list_element$f_const_mean %in% lower_outlier_ident, ]

            #Update the while loop condition
            while_loop_run <- TRUE

        } else {
            message("No lower outlier have been identified. \n")

            #Update the while loop condition
            while_loop_run <- FALSE
        }


    } else {
        message("Identifying the highest outlier on the right tail...")
        
        #An if statement to determine if there is an outlier
        if (max_z_score > crit_t) {
            
            #Assign the outlier identity and index
            upper_outlier_ident <- non_normal_list_element$f_const_mean[z_scores %in% max_z_score]
            max_outlier_index <- match(x = upper_outlier_ident, non_normal_list_element$f_const_mean)

            message("The following maximum value was found to be an outlier: ", upper_outlier_ident, "\n", "at the following row index: ", max_outlier_index, " and therefore will be removed. \n")
            
            #Assign the modified dataframe
            output_df <- non_normal_list_element[!non_normal_list_element$f_const_mean %in% upper_outlier_ident, ]

            #Update the while loop condition
            while_loop_run <- TRUE
  
        } else {
            message("No upper outlier have been identified. \n")

            #Update the while loop condition
            while_loop_run <- FALSE
        }


    }

    #Assign the output values to the output list
    output_list[[1]] <- output_df
    output_list[[2]] <- while_loop_run


    ##Return the output list
    return(output_list)

}




## Modified z-score test (non-parametric) - this test is basically the Grubb's test, just with the median and MAD instead of mean and SD
## The function identifies outliers from the input df based on the f-constants
detect_outlier_f_const_mZscore_test = function(non_normal_list_element, left_tail = TRUE) {
    ##Define function variables

    #Define the output list
    output_list <- vector(mode = "list", length = 2)
    
    #Define the f_constant vector for sorting
    ordered_f_constants <- vector(mode = "numeric", length = nrow(non_normal_list_element))
    
    #Define the median f-constant
    f_const_median <- vector(mode = "numeric", length = 1)

    #Define the median absolute deviation (MAD)
    f_const_mad <- vector(mode = "numeric", length = 1)

    #Define a vector for the z-scores
    z_scores <- vector(mode = "numeric", length = nrow(non_normal_list_element))

    #Define a vector for the minimum z-score
    min_z_score <- vector(mode = "numeric", length = 1)

    #Define a vector for the maximum z-score
    max_z_score <- vector(mode = "numeric", length = 1)

    #Define the critical z-value
    crit_t <- vector(mode = "numeric", length = 1)

    #Max outlier identity
    upper_outlier_ident <- vector(mode = "numeric", length = 1)

    #Max outlier index
    upper_outlier_index <- vector(mode = "numeric", length = 1)

    #Min outlier identity
    lower_outlier_ident <- vector(mode = "numeric", length = 1)

    #Min outlier index
    lower_outlier_index <- vector(mode = "numeric", length = 1)

    #Define the degrees of freedom
    deg_of_f <- vector(mode = "numeric", length = 1)


    ##Calculate the defined variables

    #Assign and sort the f_constants
    ordered_f_constants <- non_normal_list_element$f_const_mean[order(non_normal_list_element$f_const_mean)]

    #Calculate the median of the f-constants
    f_const_median <- median(ordered_f_constants, na.rm = TRUE)

    #Calculate the MAD
    f_const_mad <- mad(x = non_normal_list_element$f_const_mean, na.rm = TRUE)

    #Calculate the minimum z-score
    min_z_score <- (f_const_median - min(non_normal_list_element$f_const_mean)) / f_const_mad

    #Calculate all the z-scores
    for (i in seq_along(non_normal_list_element$f_const_mean)) {
        z_scores[i] <- abs((non_normal_list_element$f_const_mean[i] - f_const_median) / f_const_mad)
    
    }

    #Calculate the maximum z-score
    max_z_score <- (max(non_normal_list_element$f_const_mean) - f_const_median) / f_const_mad

    #Calculate the degrees of freedom
    deg_of_f <- nrow(non_normal_list_element) - 1

    #Calculate the critical value of t
    crit_t <- qt(1 - 0.05 / 2, deg_of_f)


    ##Compare to critical value
    
    #An if statement to control on which tail the outlier was identified
    if (left_tail == TRUE) {
        message("Identifying the lowest outlier on the left tail...")
        
        #An if statement to determine if there is an outlier
        if (min_z_score > crit_t) {
            
            #Assign the outlier identity and index
            lower_outlier_ident <- non_normal_list_element$f_const_mean[z_scores %in% min_z_score]
            lower_outlier_index <- match(x = lower_outlier_ident, non_normal_list_element$f_const_mean)

            message("The following minimum value was found to be an outlier: ", lower_outlier_ident, "at the following row index: ", lower_outlier_index, "\n")

            #Assign the detected min outlier to the output list
            output_list[[1]] <- lower_outlier_ident
            output_list[[2]] <- lower_outlier_index

            #Name the output list elements
            names(output_list) <- c("left-tail outlier", "element index")

            ##Return the output list
            return(output_list)

        } else {
            message("No lower outlier have been identified. \n")

        }


    } else {
        message("Identifying the highest outlier on the right tail...")
        
        #An if statement to determine if there is an outlier
        if (max_z_score > crit_t) {
            
            #Assign the outlier identity and index
            upper_outlier_ident <- non_normal_list_element$f_const_mean[z_scores %in% max_z_score]
            upper_outlier_index <- match(x = upper_outlier_ident, non_normal_list_element$f_const_mean)

            message("The following maximum value was found to be an outlier: ", upper_outlier_ident, "\n", "at the following row index: ", upper_outlier_index, "\n")
            
            #Assign the detected min outlier to the output list
            output_list[[1]] <- upper_outlier_ident
            output_list[[2]] <- upper_outlier_index

            #Name the output list elements
            names(output_list) <- c("right-tail outlier", "element index")

            ##Return the output list
            return(output_list)

        } else {
            message("No upper outlier have been identified. \n")

        }


    }


    

}




## A parametric outlier test function using the Grubbs test
## The function identifies adn removes outliers from the input df based on the f-constants
remove_outlier_f_const_Grubbs = function(normal_list_element, left_tail = TRUE) {
    
    
    ##Declare function variables

    #Define f-constants vector
    f_constants <- vector(mode = "numeric", length = nrow(normal_list_element))

    #Define an output dataframe
    output_df <- data.frame(matrix(nrow = nrow(normal_list_element), ncol = ncol(normal_list_element)))

    #Define the output list
    output_list <- vector(mode = "list", length = 2)

    #Grubb's test - lower outlier result
    grubbs_lower_res <- vector(mode = "list", length = 5)

    #Lower outlier identity
    lower_outlier_ident <- vector(mode = "numeric", length = 1)

    #Lower outlier index
    lower_outlier_index <- vector(mode = "numeric", length = 1)

    #Grubb's test - upper outlier result
    grubbs_upper_res <- vector(mode = "list", length = 5)

    #Upper outlier identity
    upper_outlier_ident <- vector(mode = "numeric", length = 1)

    #Upper outlier index
    upper_outlier_index <- vector(mode = "numeric", length = 1)

    #Variable to track if no outliers are found
    while_loop_run <- TRUE


    ##Calculate and assign the appropriate variables

    #Assign the f-constants to a new vector
    f_constants <- normal_list_element$f_const_mean

    #Assign the input dataframe as an output dataframe so in case of no outliers the original df will be returned
    output_df <- normal_list_element
    

    #This if statement controls if the upper or lower outlier is identified and removed
    if (left_tail == TRUE) {
        
        #Assign the value of Grubbs test on the left tail
        grubbs_lower_res <- outliers::grubbs.test(f_constants, opposite = TRUE)

        #This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_lower_res$p.value < 0.05) {

            #Assign outlier identity and index
            lower_outlier_ident <- as.numeric(stringr::str_extract(grubbs_lower_res$alternative, "[0-9]+\\.[0-9]+"))
            lower_outlier_index <- match(x = lower_outlier_ident, f_constants)

            #An if statement to make sure that the output is only returned if there was an outlier    
            if (length(lower_outlier_index) > 0) {
                message("The following lower outliers will be removed: ", lower_outlier_ident, "\n")
                
                #Assign the modified dataframe
                output_df <- normal_list_element[!f_constants %in% lower_outlier_ident, ]

                #Update the while loop condition
                while_loop_run <- TRUE
                
            }


        } else {
            message("No significant outlier found on the left tail. \n")
            
            #Update the while loop condition
            while_loop_run <- FALSE
            

        }


    } else {
        
        #Assign the value of Grubbs test on the right tail
        grubbs_upper_res <- outliers::grubbs.test(f_constants)

        #This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_upper_res$p.value < 0.05) {
            
            #Assign outlier identity and index
            upper_outlier_ident <- as.numeric(stringr::str_extract(grubbs_upper_res$alternative, "[0-9]+\\.[0-9]+"))
            upper_outlier_index <- match(x = upper_outlier_ident, f_constants)
                
            #An if statement to make sure that the output is only returned if there was an outlier
            if (length(upper_outlier_index) > 0) {
                message("The following upper outliers will be removed: ", upper_outlier_ident, "\n")
                
                #Assign the modified dataframe
                output_df <- normal_list_element[!f_constants %in% upper_outlier_ident, ]

                #Update the while loop condition
                while_loop_run <- TRUE
                
            }


        } else {
            message("No significant outlier found on the right tail. \n")
            
            #Update the while loop condition
            while_loop_run <- FALSE
            
        }
    }


    #Assign the output values to the output list
    output_list[[1]] <- output_df
    output_list[[2]] <- while_loop_run


    ##Return the output list
    return(list(output_df, while_loop_run))
    
}




## A parametric outlier test function using the Grubbs test
## The function identifies adn removes outliers from the input df based on the f-constants
detect_outlier_f_const_Grubbs = function(normal_list_element, left_tail = TRUE) {  
    ##Declare function variables

    #Define f-constants vector
    f_constants <- vector(mode = "numeric", length = nrow(normal_list_element))

    #Define the output list
    output_list <- vector(mode = "list", length = 2)

    #Grubb's test - lower outlier result
    grubbs_lower_res <- vector(mode = "list", length = 5)

    #Lower outlier identity
    lower_outlier_ident <- vector(mode = "numeric", length = 1)

    #Lower outlier index
    lower_outlier_index <- vector(mode = "numeric", length = 1)

    #Grubb's test - upper outlier result
    grubbs_upper_res <- vector(mode = "list", length = 5)

    #Upper outlier identity
    upper_outlier_ident <- vector(mode = "numeric", length = 1)

    #Upper outlier index
    upper_outlier_index <- vector(mode = "numeric", length = 1)


    ##Calculate and assign the appropriate variables

    #Assign the f-constants to a new vector
    f_constants <- normal_list_element$f_const_mean

    #This if statement controls if the upper or lower outlier is identified and removed
    if (left_tail == TRUE) {
        
        #Assign the value of Grubbs test on the left tail
        grubbs_lower_res <- outliers::grubbs.test(f_constants, opposite = TRUE)

        #This if statement controls that if there is an outlier it should be identified and removed
        if (grubbs_lower_res$p.value < 0.05) {

            #Assign outlier identity and index
            lower_outlier_ident <- as.numeric(stringr::str_extract(grubbs_lower_res$alternative, "[0-9]+\\.[0-9]+"))
            lower_outlier_index <- match(x = lower_outlier_ident, f_constants)

            #An if statement to make sure that the output is only returned if there was an outlier    
            if (length(lower_outlier_index) > 0) {
                message("The following minimum value was found to be an outlier: ", lower_outlier_ident, "at the following row index: ", lower_outlier_index)
                
                #Assign the detected min outlier to the output list
                output_list[[1]] <- lower_outlier_ident
                output_list[[2]] <- lower_outlier_index

                #Name the output list elements
                names(output_list) <- c("left-tail outlier", "element index")

                #Return the output list
                return(output_list)
                
            }


        } else {
            message("No significant outlier found on the left tail.")

        }

    } else {
        
        #Assign the value of Grubbs test on the right tail
        grubbs_upper_res <- outliers::grubbs.test(f_constants)

        #This if statement controls that if there is an outlier it should be identified
        if (grubbs_upper_res$p.value < 0.05) {
            
            # Assign outlier identity and index
            upper_outlier_ident <- as.numeric(stringr::str_extract(grubbs_upper_res$alternative, "[0-9]+\\.[0-9]+"))
            upper_outlier_index <- match(x = upper_outlier_ident, f_constants)
                
            #An if statement to make sure that the output is only returned if there was an outlier
            if (length(upper_outlier_index) > 0) {
                message("The following maximum value was found to be an outlier: ", upper_outlier_ident, "\n", "at the following row index: ", upper_outlier_index)
                
                #Assign the detected min outlier to the output list
                output_list[[1]] <- upper_outlier_ident
                output_list[[2]] <- upper_outlier_index

                #Name the output list elements
                names(output_list) <- c("right-tail outlier", "element index")

                #Return the output list
                return(output_list)
                
            }

        } else {
            message("No significant outlier found on the right tail.")
            
        }
    }

    
    
}




## This a wrapper function which binds the 3 outlier detection and removal functions together and runs them in a  while loop until there are no more
## outliers
## IMPORTANT NOTE: this function seems to work, but is a bit cobbled together as I'm not very familiar with while loops. Therefore please if you have more
## experience with them check and debug this bad boy! :)
outlier_cleaner = function(calculate_f_constants_output_list, is_data_normal_output_list, nonparam_test = "numeric_outlier_test") {
    for (index in seq_along(is_data_normal_output_list)) {
        if (is_data_normal_output_list[[index]]$p.value < 0.05) {
            message("\nIt seems your f-constants show a non-normal distribution. The chosen non-parametric test will be used for outlier detection and removal.")

            if (nonparam_test == "numeric_outlier_test") {
                #Print the name of the current list element for better orientation
                message(names(calculate_f_constants_output_list)[[index]])

                #Run the appropriate outlier cleaner function and assign the result to the appropriate output list element
                calculate_f_constants_output_list[[index]] <- remove_outlier_f_const_NOTest(calculate_f_constants_output_list[[index]])
                
            } else if (nonparam_test == "mZscore_test") {
                #Print the name of the current list element for better orientation
                message(names(calculate_f_constants_output_list)[[index]])

                #Set the eval condition for the while loop
                while_loop_run <- TRUE

                #Initiate a while loop which will loop through the dataframe while outliers on the left tail are present
                while (while_loop_run == TRUE) {
                    #Run the appropriate outlier cleaner function and assign the result to a new list which will hold the two elements
                    result <- remove_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[index]], left_tail = TRUE)

                    #Assign the appropriate outlier cleaner function's result to the appropriate output list element
                    calculate_f_constants_output_list[[index]] <- result[[1]]

                    #Update while_loop_run's eval parameter based on outliers_found
                    while_loop_run <- result[[2]]  
                    
                }

                #Set the eval condition for the while loop
                while_loop_run <- TRUE

                #Initiate a while loop which will loop through the dataframe while outliers on the right tail are present
                while (while_loop_run == TRUE) {
                    #Run the appropriate outlier cleaner function and assign the result to a new list which will hold the two elements
                    result <- remove_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[index]], left_tail = FALSE)

                    #Assign the appropriate outlier cleaner function's result to the appropriate output list element
                    calculate_f_constants_output_list[[index]] <- result[[1]]

                    #Update while_loop_run's eval parameter based on outliers_found
                    while_loop_run <- result[[2]]
                    
                }

            }


        } else {
            message("\nIt seems your f-constants show a normal distribution. A parametric Grubbs test will be used for outlier detection and removal.")

            #Print the name of the current list element for better orientation
            message(names(calculate_f_constants_output_list)[[index]])

            #Set the eval condition for the while loop
            while_loop_run <- TRUE
            
            #Initiate a while loop which will loop through the dataframe while outliers on the left tail are present
            while (while_loop_run == TRUE) {
                #Run the appropriate outlier cleaner function and assign the result to a new list which will hold the two elements
                result <- remove_outlier_f_const_Grubbs(calculate_f_constants_output_list[[index]], left_tail = TRUE)

                #Assign the appropriate outlier cleaner function's result to the appropriate output list element
                calculate_f_constants_output_list[[index]] <- result[[1]]

                #Update while_loop_run's eval parameter based on outliers_found
                while_loop_run <- result[[2]]
                    
            }

            #Set the eval condition for the while loop
            while_loop_run <- TRUE

            #Initiate a while loop which will loop through the dataframe while outliers on the right tail are present
            while (while_loop_run == TRUE) {
                #Run the appropriate outlier cleaner function and assign the result to a new list which will hold the two elements 
                result <- remove_outlier_f_const_Grubbs(calculate_f_constants_output_list[[index]], left_tail = FALSE)

                #Assign the appropriate outlier cleaner function's result to the appropriate output list element
                calculate_f_constants_output_list[[index]] <- result[[1]]

                #Update while_loop_run's eval parameter based on outliers_found
                while_loop_run <- result[[2]]
                    
            }

        }
    }


    ## Return the outlier free output list
    return(calculate_f_constants_output_list)
}




## This a wrapper function which binds the 3 outlier detection functions together (please note that the detect function returns nothing, only messages!)
outlier_detector = function(calculate_f_constants_output_list, is_data_normal_output_list, nonparam_test = "numeric_outlier_test") {
    ##Run the appropriate detector functions

    #This construct checks if the f-constants from the elements of the calculate f-constants output list are normally distributed
    #NOTE: this is possible because the element order is preserved between the calc.f-const list and the is_data_normal list 
    for (index in seq_along(is_data_normal_output_list)) {
        
        #This if statement checks if the f-factor distribution is normal
        if (is_data_normal_output_list[[index]]$p.value < 0.05) {
            message("\nIt seems your f-constants show a non-normal distribution. The chosen non-parametric test will be used for outlier detection.")
        
            #This if statement is responsible for choosing a non-parametric test
            if (nonparam_test == "numeric_outlier_test") {
                
                #Print the name of the current list element for better orientation
                message(names(calculate_f_constants_output_list)[[index]])

                #Run the appropriate detection test
                detect_outlier_f_const_NOTest(calculate_f_constants_output_list[[index]])
  
            } else if (nonparam_test == "mZscore_test") {

                
                #Print the name of the current list element for better orientation
                message(names(calculate_f_constants_output_list)[[index]])

                #Run the appropriate detection test
                detect_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[index]], left_tail = TRUE)

                #Run the appropriate detection test
                detect_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[index]], left_tail = FALSE)

            }
        
        } else {
            message("\nIt seems your f-constants show a normal distribution. A parametric Grubbs test will be used for outlier detection.")

            
            #Print the name of the current list element for better orientation
            message(names(calculate_f_constants_output_list)[[index]])

            #Run the appropriate detection test
            detect_outlier_f_const_Grubbs(calculate_f_constants_output_list[[index]], left_tail = TRUE)

            #Run the appropriate detection test
            detect_outlier_f_const_Grubbs(calculate_f_constants_output_list[[index]], left_tail = FALSE)
                    
        }

    }   

}


################################################# Section end #################################################





################################################    MARK: calculate the mean f-constants    #################################################
                                                #       and estimate tumor volumes          #




## Thi function will calculate the mean  f-constants
calc_mean_f = function(calculate_f_constants_output_list) {
    ##Declare the function variables

    #List to return
    return_list <- vector(mode = "list", length = length(calculate_f_constants_output_list))


    ##Assign function variables
    input_list <- calculate_f_constants_output_list


    ##Calculations
    for (index in seq_along(input_list)) {
        return_list[[index]] <- mean(input_list[[index]]$f_const_mean)
    }

    #Name the return list elements
    names(return_list) <- names(input_list)

    ##Return the output list
    return(return_list)

}




## This function will estimate the tumor volumes based on the mean f-constant
estimate_tumor_volume = function(input_measurement_list, mean_f_values_list) {
    ##Declare the function variables

    #List to return
    return_list <- vector(mode = "list", length = length(input_measurement_list))


    ##Assign function variables
    input_list <- input_measurement_list


    ##Calculations
    
    #This for loop will walk along the input list to access each element
    #NOTE: the length of the input_df and the mean_f_values_list is the same, so indexes can be used for both
    for (index in seq_along(input_list)) {
        #Initialize a temp_df using the elements of the input list
        temp_df <- input_list[[index]]

        #Initialize a vector containing the sampling dates for naming the estimated tumor volumes
        sampling_dates <- unique(unlist(stringr::str_extract_all(string = colnames(temp_df), pattern = "_[0-9\\.]+")))

        #Initialize a temporary df which will contain the LxW columns for the volume estimation
        temp_lxw_df <- temp_df[, grep(pattern = "LxW", x = colnames(temp_df), ignore.case = TRUE)]

        #Initialize a temporary df which will contain the estimated tumor volumes
        temp_vol_df <- data.frame(matrix(nrow = nrow(temp_lxw_df), ncol = ncol(temp_lxw_df)))

        #This nested for loop will progress through the LxW columns and rows of each individual list element (dataframe) to estimate the tumor volumes using the 
        #appropriate formula
        for (col in seq_len(ncol(temp_lxw_df))) {
            for (row in seq_len(nrow(temp_lxw_df))) {
            
                #Estimate and assign the tumor volumes to each individual sample in the temp_df
                temp_vol_df[row, col] <- (pi / 6) * mean_f_values_list[[index]] * (temp_lxw_df[row, col])^(3 / 2)
            }
        }

        #Name the new columns containing the estimated tumor volumes
        colnames(temp_vol_df) <- paste0("Estim_volume", sampling_dates)
        
        #Bind the estimated volumes to the parent dataframe
        temp_df <- cbind(temp_df, temp_vol_df)

        #Assign the temp_df to the return list for each input list element
        return_list[[index]] <- temp_df

    }

    #Name the return_list elements according to the input list
    names(return_list) <- names(input_measurement_list)


    ## Return the modified dataframes
    return(return_list)

}


################################################# Section end #################################################





################################################    MARK: correct the estimated volumes    #################################################
                                                #       and test the goodness of fit       #



###MARK: Cntinue from here
## Correct the estimated tumor volumes
## NOTE: this is a complex function. The idea is that it tests what percentage of tumor volume deviation comes from the 0.5, 1, 1.5, 2 
## standard deviation distance of sample f-constants from the calculated mean f-constant. Then based on which SD bracket
##f-constants fall into it corrects the estimated tumor volume by the appropriate percentage depending on if the value
## is under or over estimated (this depends on which direction the f-constant is deviating from the mean)
tumor_vol_correction = function(estimated_tumor_volume_list, mean_f_values_list) {
    ##Declare static function variables

    #Initialize the return/output list
    return_list <- vector(mode = "list", length = length(estimated_tumor_volume_list))

    correction_factor_list <- vector(mode = "list", length = length(estimated_tumor_volume_list))


    ##Assign the appropriate variables

    temp_df <- data.frame()

    #Assign the input list
    input_list <- estimated_tumor_volume_list


    ##Calculate the appropriate variables

    for (index in seq_along(input_list)) {
        ##Declare dynamic function variables

        #Define a vector for the standard deviation bins
        std_bin <- vector(mode = "numeric", length = 5)

        #Define a vector containing the standard deviations of 0.5-2 SD below mean
        below_mean_std_vals <- vector(mode = "numeric", length = 5)
        names(below_mean_std_vals) <- as.character(seq(from = 0, to = 2, by = 0.5))

        #Define a vector containing the standard deviations of 0.5-2 SD above mean
        above_mean_std_vals <- vector(mode = "numeric", length = 5)
        names(above_mean_std_vals) <- as.character(seq(from = 0, to = 2, by = 0.5))

        #Define a vector containing the tumor volume deviation (in percentage) of 0.5-2 SD below mean
        below_mean_volumes <- data.frame(matrix(nrow = 5, ncol = 2))
        rownames(below_mean_volumes) <- as.character(seq(from = 0, to = 2, by = 0.5))
        colnames(below_mean_volumes) <- paste0("vol_deviation", uCT_sampling_dates)

        #Define a vector containing the tumor volume deviation (in percentage) of 0.5-2 SD above mean
        above_mean_volumes <- data.frame(matrix(nrow = 5, ncol = 2))
        rownames(above_mean_volumes) <- as.character(seq(from = 0, to = 2, by = 0.5))
        colnames(below_mean_volumes) <- paste0("vol_deviation", uCT_sampling_dates)

        #Define a dynamic filter_df
        filter_list <- vector(mode = "list", length = 4)

        #Define the standard deviation variable
        temp_std <- vector(mode = "numeric", length = 1)
        

        ##Assign the dynamic function variables

        #Assign the temp_df 
        temp_df <- input_list[[index]]

        estim_vols <- temp_df[, grep(pattern = "estim", x = colnames(temp_df), ignore.case = TRUE)]

        uCT_vols <- temp_df[, grep(pattern = "uct", x = colnames(temp_df), ignore.case = TRUE)]

        #Assign the standard deviation variable
        temp_std <- sd(temp_df$f_const_mean)

        #Assign the standard deviation fractions
        std_bin <- seq(from = 0, to = 2, by = 0.5)

        #Assign the standard deviations below mean (we start with 0 - the mean - and go down by half a sd compared to the mean)
        below_mean_std_vals <- mean_f_values_list[[index]] - (temp_std * std_bin)
        
        #Assign the standard deviations above mean (we start with 0 - the mean - and go up by half a sd compared to the mean)
        above_mean_std_vals <- mean_f_values_list[[index]] + (temp_std * std_bin)

        #Assign the tumor volume deviation below mean (in this case the tumor volumes are
        #over-estimated)
        for (c in seq_len(ncol(estim_vols))) {
            for (ind in seq_along(below_mean_std_vals)) {
                filter_list[[ind]] <- dplyr::filter(temp_df, f_const_mean <= below_mean_std_vals[ind] & f_const_mean >= below_mean_std_vals[ind + 1])
                
                below_mean_volumes[ind, c] <- mean(filter_list[[ind]][, colnames(estim_vols)[c]], na.rm = TRUE) / mean(filter_list[[ind]][, colnames(uCT_vols)[c]], na.rm = TRUE)
            
            }

        }
        
        #Assign the tumor volume deviation above mean (in this case the tumor volumes are
        #under-estimated)
        for (c in seq_along(estim_vols)) {
            for (ind in seq_along(above_mean_std_vals)) {
                filter_list[[ind]] <- dplyr::filter(temp_df, f_const_mean > above_mean_std_vals[ind] & f_const_mean < above_mean_std_vals[ind + 1])
            
                above_mean_volumes[ind, c] <- mean(filter_list[[ind]][, colnames(estim_vols)[c]], na.rm = TRUE) / mean(filter_list[[ind]][, colnames(uCT_vols)[c]], na.rm = TRUE)
            
            }

        }
        
        #corrected_volumes <- vector(mode = "numeric", length = nrow(temp_df))
        corrected_volumes <- data.frame(matrix(nrow = nrow(temp_df), ncol = 2))
        correction_factors <- data.frame(matrix(nrow = nrow(temp_df), ncol = 2))

        #Perform the value correction for values which f-constants are below the mean f-constant (in this case the tumor volumes are
        #over-estimated). The outer for loop traverses the columns of the split off estimated_volumes dfs
        for (c in seq_len(ncol(estim_vols))) {
            
            #The inner for loop traverses the columns of the temp_df (note: this is interchangeable with the estim_vols df as they have the same number of rows)
            for (r in seq_len(nrow(temp_df))) {
 
                #The innermost for-loop traverses the standard deviation values of f-constants below the mean f-constant (the SD here is calculated for the 
                #f-constants and brackets of 0(mean)-0.5-1-1.5-2 SD are used)
                for (e in seq_len(length(below_mean_std_vals) - 1)) {

                    #An if statement to control which correction factor (this is volume parentage of the estimated volume compared to the measured volume)
                    #should be used, based on which f-constant bracket the sample falls into
                    if (temp_df$f_const_mean[r] <= below_mean_std_vals[e] && temp_df$f_const_mean[r] >= below_mean_std_vals[e + 1]) {

                        #Perform the correction by dividing the estimated tumor volume with the correction factor (the percentage of over or under estimation
                        #as a function of the f-constant deviation) and assign the corrected values into a new dataframe
                        corrected_volumes[r, c] <- temp_df[r, colnames(estim_vols)[c]] / below_mean_volumes[e, c]
                        correction_factors[r, c] <- below_mean_volumes[e, c]

                    }
                
                }

            }
            
        }

        #Perform the value correction for values which f-constants are above the mean f-constant (in this case the tumor volumes are
        #under-estimated). The outer for loop traverses the columns of the split off estimated_volumes dfs
        for (c in seq_len(ncol(estim_vols))) {
            
            #The inner for loop traverses the columns of the temp_df (note: this is interchangeable with the estim_vols df as they have the same number of rows)
            for (r in seq_len(nrow(temp_df))) {

                #The innermost for-loop traverses the standard deviation values of f-constants below the mean f-constant (the SD here is calculated for the 
                #f-constants and brackets of 0(mean)-0.5-1-1.5-2 SD are used)
                for (e in seq_len(length(above_mean_std_vals) - 1)) {

                    #An if statement to control which correction factor (this is volume parentage of the estimated volume compared to the measured volume)
                    #should be used, based on which f-constant bracket the sample falls into
                    if (temp_df$f_const_mean[r] >= above_mean_std_vals[e] && temp_df$f_const_mean[r] <= above_mean_std_vals[e + 1]) {

                        #Perform the correction by dividing the estimated tumor volume with the correction factor (the percentage of over or under estimation
                        #as a function of the f-constant deviation) and assign the corrected values into a new daraframe
                        corrected_volumes[r, c] <- temp_df[r, colnames(estim_vols)[c]] / above_mean_volumes[e, c]
                        correction_factors[r, c] <- above_mean_volumes[e, c]

                    }
                
                }

            }
             
        }

        #Name the columns of the correction factor dataframe
        colnames(correction_factors) <- paste0("Corr_factor", uCT_sampling_dates)

        #Assign the correction factor dataframes to the correction factor return list
        correction_factor_list[[index]] <- correction_factors
        

        #Name the columns of the corrected volumes dataframe
        colnames(corrected_volumes) <- paste0("Corr_volumes", uCT_sampling_dates)

        #Bind the corrected volumes dataframe to the temp_df
        temp_df <- cbind(temp_df, corrected_volumes)

        #Assign the updated temp_df to the return list
        return_list[[index]] <- temp_df

    }

    #Name the correction_factor_list elements according to the input list
    names(correction_factor_list) <- names(estimated_tumor_volume_list)

    #Assign the correction factor list to the global environment with the sam name (an alternative form of return)
    correction_factor_list <<- correction_factor_list

    #Name the return_list elements according to the input list
    names(return_list) <- names(estimated_tumor_volume_list)

    
    ## Return the output list
    return(return_list)

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
