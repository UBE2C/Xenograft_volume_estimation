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
read_uCT_data = function(data_path) {
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
read_caliper_data = function(data_path) {
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


#MARK: function in work continue from here!!!
fit_caliper_measurements = function(read_caliper_data_output_list, clean_uCT_data_output_list) {
    ##Define the variables used by this function
    
    #Initialize a list onto which the row trimmed down caliper measurements will be added
    trimmed_caliper_list <- vector(mode = "list")

    #Initialize a list onto which the column trimmed down caliper measurements will be added
    col_trimmed_caliper_list <- list()

    #Initialize a list onto which the grouped trimmed down caliper measurements will be added
    grouped_trimmed_caliper_list <- list()
    
    #Initialize a temporary dataframe on which the processing can take place
    temp_df <- data.frame()

    #Initialize temporary date lists
    temp_date_list <- list()

    #Initialize a list onto which the unique dates of the corresponding caliper measurements will be added
    corresponding_dates <- vector(mode = "list", length = length(clean_uCT_data_output_list))

    #Define the treatment group_IDs shared by both uCT and caliper_measurement files
    group_IDs <- vector(mode = "list", length = length(clean_uCT_data_output_list))

    for (i in seq_along(clean_uCT_data_output_list)) {
        group_IDs[[i]] <- unique(clean_uCT_data_output_list[[i]]$Treatment_group_ID)
    }


    ##Start the processing steps

    #First, trim the caliper measurements by sample/mouse ID to match the uCT measurements

    #A for loop which trims down the caliper measurements dataframes to the corresponding uCT measurements Mouse_IDs
    #The outer for loop traverses the caliper measurements list of dataframes
    for (i in seq_along(read_caliper_data_output_list)) {

        #The inner for loop traverses the list of uCT measurements dataframes 
        for (e in seq_along(clean_uCT_data_output_list)) {

            #This line filters the caliper measurements by the mouse_ID found in the corresponding uCT dataframes
            temp_df <- dplyr::filter(read_caliper_data_output_list[[i]], read_caliper_data_output_list[[i]]$Mouse_ID %in% clean_uCT_data_output_list[[e]]$Mouse_ID)

            #This if statement makes sure that only non 0 row dataframes gets assigned adn assigns them to a list
            if (nrow(temp_df) > 0) {

                #NOTE: the use of value = list() is necessary as This is necessary because append() expects the new value to be a list itself.
                #If temp_df were not a list, it would be automatically coerced to one and then the formatting is all wrong
                trimmed_caliper_list <- append(x = trimmed_caliper_list, values = list(temp_df))
            }
            
        }

        
    }

    #A for loop which traverses the clean uCT measurement list to assign the corresponding dates to the proper list
    for (i in seq_along(clean_uCT_data_output_list)) {
            
            #An if statements which controls which date columns will be used for sample matching
            if (is.element(el = "Corresponding_caliper_measurement_dates", set = colnames(x = clean_uCT_data_output_list[[i]]))) {
                
                #Assign the corresponding dates to the proper list
                temp_date_list[[i]] <- unique(clean_uCT_data_output_list[[i]]$Corresponding_caliper_measurement_dates)
            } else {
                warning("The column: 'Corresponding_caliper_measurement_dates' was not found in the", i, "element of the uCT data list.",
                        "\n  Assuming that all caliper measurements correspond with the uCT measurements, the standard miCT_dates column will be used.")
                
                #Assign the dates to the proper list
                temp_date_list[[i]] <- unique(clean_uCT_data_output_list[[i]]$miCT_dates)
                
            }
    }

    #This for loop splits up the corresponding dates (separated by "_") and uses unlist to store them as vectors in a list
    for (i in seq_along(temp_date_list)) {
        corresponding_dates[[i]] <- unlist(stringr::str_split(temp_date_list[[i]], pattern = "_"))
    }

    #This for loop converts char "NA"s to real NAs, removes these from the vectors and sorts out the unique dates
    for (i in seq_along(corresponding_dates)) {
        corresponding_dates[[i]] <- dplyr::na_if(x = corresponding_dates[[i]], y = "NA")
        corresponding_dates[[i]] <- corresponding_dates[[i]][!is.na(corresponding_dates[[i]])]
        corresponding_dates[[i]] <- unique(corresponding_dates[[i]])
    }

    #Assign the corresponding_dates to the parent main() environment
    assign("corresponding_dates_main", corresponding_dates, envir = parent.frame())

    #Select the corresponding LxW measurements from the trimmed caliper measurements list to the uCT list
    for (i in seq_along(trimmed_caliper_list)) {
        
        #Add the first 3 ID columns to the empty list elements so the LxW columns can be bound to them next
        col_trimmed_caliper_list[[i]] <- trimmed_caliper_list[[i]][, c("Mouse_ID", "Treatment_group_ID", "Treatment")]
        
        #This for loop binds the LxW columns which correspond to the uCT measurements
        for (e in seq_along(corresponding_dates)) {
            col_trimmed_caliper_list[[i]] <- cbind(col_trimmed_caliper_list[[i]], dplyr::select(.data = trimmed_caliper_list[[i]], matches(corresponding_dates[[e]])))
        }
    }
    
    #Remove duplicated columns from the column trimmed list's dataframes
    for (i in seq_along(col_trimmed_caliper_list)) {
       col_trimmed_caliper_list[[i]] <- col_trimmed_caliper_list[[i]][, unique(colnames(col_trimmed_caliper_list[[i]]))] 
    }

    #Initialize a dynamic list onto which the grouping variables will be added
    #grouping_list <- vector(mode = "list", length = 2)    

    #Group the selected columns according to uCT groups
    #for (i in seq_along(corresponding_dates)) {
    #    grouping_list[[i]] <- dplyr::select(.data = temp_df, matches(corresponding_dates[[i]]))
    #}
    

    
    
    return(col_trimmed_caliper_list)

    
}


#A for loop which traverses the clean uCT measurement list to assign the corresponding dates to the proper list
    for (i in seq_along(clean_uCT_data_output_list)) {
            
            #An if statements which controls which date columns will be used for sample matching
            if (is.element(el = "Corresponding_caliper_measurement_dates", set = colnames(x = clean_uCT_data_output_list[[i]]))) {
                
                #Assign the corresponding dates to the proper list
                temp_date_list[[i]] <- unique(clean_uCT_data_output_list[[i]]$Corresponding_caliper_measurement_dates)
            } else {
                warning("The column: 'Corresponding_caliper_measurement_dates' was not found in the", i, "element of the uCT data list.",
                        "\n  Assuming that all caliper measurements correspond with the uCT measurements, the standard miCT_dates column will be used.")
                
                #Assign the dates to the proper list
                temp_date_list[[i]] <- unique(clean_uCT_data_output_list[[i]]$miCT_dates)
                
            }
    }

    #This for loop splits up the corresponding dates (separated by "_") and uses unlist to store them as vectors in a list
    for (i in seq_along(temp_date_list)) {
        corresponding_dates[[i]] <- unlist(stringr::str_split(temp_date_list[[i]], pattern = "_"))
    }

    #This for loop converts char "NA"s to real NAs, removes these from the vectors and sorts out the unique dates
    for (i in seq_along(corresponding_dates)) {
        corresponding_dates[[i]] <- dplyr::na_if(x = corresponding_dates[[i]], y = "NA")
        corresponding_dates[[i]] <- corresponding_dates[[i]][!is.na(corresponding_dates[[i]])]
        corresponding_dates[[i]] <- unique(corresponding_dates[[i]])
    }

    #Assign the corresponding_dates to the parent main() environment
    assign("corresponding_dates_main", corresponding_dates, envir = parent.frame())

    #Select the corresponding LxW measurements from the caliper measurements list to the uCT list
    for (i in seq_along(read_caliper_data_output_list)) {
        for (e in seq_along(corresponding_dates)) {
            trimmed_caliper_list <- append(trimmed_caliper_list, dplyr::select(.data = read_caliper_data_output_list[[i]], matches(corresponding_dates[[e]])))
        }
    }

    #Bind the trimmed caliper measurements columns bach into a dataframe for easier handling
    temp_df <- as.data.frame(do.call(cbind, trimmed_caliper_list))
    
    #Remove duplicated columns from the temp dataframe
    temp_df <- temp_df[, unique(colnames(temp_df))]

    #Initialize a dynamic list onto which the grouping variables will be added
    grouping_list <- vector(mode = "list", length = 2)    

    #Group the selected columns according to uCT groups
    for (i in seq_along(corresponding_dates)) {
        grouping_list[[i]] <- dplyr::select(.data = temp_df, matches(corresponding_dates[[i]]))
    }
    





for (i in seq_along(clean_uCT_list)) {
        for (e in seq_along(calip_data)) {
            print(calip_data[[e]][clean_uCT_list[[i]]$Mouse_ID %in% calip_data[[e]]$Mouse_ID, ])
        }
    }

lst <- vector(mode = "list", length = 4)
for (i in seq_along(calip_data)) {
    for (e in seq_along(clean_uCT_list)){
        print(calip_data[[i]][calip_data[[i]]$Mouse_ID %in% clean_uCT_list[[e]]$Mouse_ID, ])   
    }
}

lst <- vector(mode = "list", length = length(calip_data))
for (i in seq_along(calip_data)) {
    for (e in seq_along(clean_uCT_list)){
        # Only assign to lst[[i]] if it's the correct index
        if (i == e) {
            lst[[i]] <- calip_data[[i]][calip_data[[i]]$Mouse_ID %in% clean_uCT_list[[e]]$Mouse_ID, ]
        }
    }
}


result_list <- list()

for (i in seq_along(calip_data)) {
    result_list[[i]] <- list()  # Initialize the sublist for each iteration
    for (e in seq_along(clean_uCT_list)){
        result_list[[i]][[e]] <- calip_data[[i]][calip_data[[i]]$Mouse_ID %in% clean_uCT_list[[e]]$Mouse_ID, ]
    }
}

dplyr::filter(calip_data[[1]], calip_data[[1]]$Mouse_ID %in% clean_uCT_list[[1]]$Mouse_ID)






calip_data[[1]][calip_data[[1]]$Mouse_ID %in% clean_uCT_list[[1]]$Mouse_ID, ]
calip_data[[1]][calip_data[[1]]$Mouse_ID %in% clean_uCT_list[[2]]$Mouse_ID, ]
calip_data[[2]][calip_data[[2]]$Mouse_ID %in% clean_uCT_list[[1]]$Mouse_ID, ]
calip_data[[2]][calip_data[[2]]$Mouse_ID %in% clean_uCT_list[[2]]$Mouse_ID, ]


calip_data[[1]][calip_data[[1]]$Mouse_ID %in% clean_uCT_list[[1]]$Mouse_ID, ]
calip_data[[2]][calip_data[[2]]$Mouse_ID %in% clean_uCT_list[[1]]$Mouse_ID, ]



















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


    ##Return the bound_df_list
    return(unified_df_list)

}


################################################# Section end #################################################





################################################    MARK: f-constant    #################################################
                                                #  calculation and eval #




## Calculate the f-constants to each uCT-measurement-caliper measurement group


# NOTE: original formula V = (pi/6) * f * (l * w)^(3/2)
# NOTE: formula solved to f f = V / (pi/6) * (l * w)^(3/2)


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


    ##Return the output list containing the f-constants
    return(measurements_with_f_constants_list)
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
        shapiro_results[[i]] <- shapiro.test(calculate_f_constants_output_list[[i]]$f_constants)
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
    ordered_constants <- list_element_to_modify$f_constants[order(list_element_to_modify$f_constants)]


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
        list_element_to_modify <- list_element_to_modify[!list_element_to_modify$f_constants %in% upper_outlier, ]
    }

    #An if statement to check if there were any lower outliers and if yes to remove them
    if (length(lower_outlier) == 0) {
        message("No outlier found on the left (lower) tail for input list element. \n")
    } else {
        message("The following elements are outliers on the left (lower) tail, and will be removed: \n", lower_outlier, "\n")

        #Remove the lower outliers
        list_element_to_modify <- list_element_to_modify[!list_element_to_modify$f_constants %in% lower_outlier, ]
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
    ordered_constants <- non_normal_list_element$f_constants[order(non_normal_list_element$f_constants)]


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
            
    }

    #An if statement to check if there were any lower outliers and if yes to remove them
    if (length(lower_outlier) == 0) {
        message("No outlier found on the left (lower) tail for input list element. \n")

    } else {
        message("The following elements are outliers on the left (lower) tail: \n", lower_outlier, "\n")

    }



    #Assign the outlier f-constants to the output list
    outliers_list[[1]] <- upper_outlier
    outliers_list[[2]] <- lower_outlier


    ##Return the output list
    return(outliers_list)
    
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
    ordered_f_constants <- non_normal_list_element$f_constants[order(non_normal_list_element$f_constants)]

    #Assign the input dataframe as an output dataframe so in case of no outliers the original df will be returned
    output_df <- non_normal_list_element

    #Calculate the median of the f-constants
    f_const_median <- median(ordered_f_constants, na.rm = TRUE)

    #Calculate the MAD
    f_const_mad <- mad(x = non_normal_list_element$f_constants, na.rm = TRUE)

    #Calculate the minimum z-score
    min_z_score <- (f_const_median - min(non_normal_list_element$f_constants)) / f_const_mad

    #Calculate all the z-scores
    for (i in seq_along(non_normal_list_element$f_constants)) {
        z_scores[i] <- abs((non_normal_list_element$f_constants[i] - f_const_median) / f_const_mad)
    
    }

    #Calculate the maximum z-score
    max_z_score <- (max(non_normal_list_element$f_constants) - f_const_median) / f_const_mad

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
            lower_outlier_ident <- non_normal_list_element$f_constants[z_scores %in% min_z_score]
            min_outlier_index <- match(x = lower_outlier_ident, non_normal_list_element$f_constants)

            message("The following minimum value was found to be an outlier: ", lower_outlier_ident, "at the following row index: ", min_outlier_index, " and therefore will be removed. \n")
            
            #Assign the modified dataframe
            output_df <- non_normal_list_element[!non_normal_list_element$f_constants %in% lower_outlier_ident, ]

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
            upper_outlier_ident <- non_normal_list_element$f_constants[z_scores %in% max_z_score]
            max_outlier_index <- match(x = upper_outlier_ident, non_normal_list_element$f_constants)

            message("The following maximum value was found to be an outlier: ", upper_outlier_ident, "\n", "at the following row index: ", max_outlier_index, " and therefore will be removed. \n")
            
            #Assign the modified dataframe
            output_df <- non_normal_list_element[!non_normal_list_element$f_constants %in% upper_outlier_ident, ]

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
    ordered_f_constants <- non_normal_list_element$f_constants[order(non_normal_list_element$f_constants)]

    #Calculate the median of the f-constants
    f_const_median <- median(ordered_f_constants, na.rm = TRUE)

    #Calculate the MAD
    f_const_mad <- mad(x = non_normal_list_element$f_constants, na.rm = TRUE)

    #Calculate the minimum z-score
    min_z_score <- (f_const_median - min(non_normal_list_element$f_constants)) / f_const_mad

    #Calculate all the z-scores
    for (i in seq_along(non_normal_list_element$f_constants)) {
        z_scores[i] <- abs((non_normal_list_element$f_constants[i] - f_const_median) / f_const_mad)
    
    }

    #Calculate the maximum z-score
    max_z_score <- (max(non_normal_list_element$f_constants) - f_const_median) / f_const_mad

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
            lower_outlier_ident <- non_normal_list_element$f_constants[z_scores %in% min_z_score]
            lower_outlier_index <- match(x = lower_outlier_ident, non_normal_list_element$f_constants)

            message("The following minimum value was found to be an outlier: ", lower_outlier_ident, "at the following row index: ", lower_outlier_index, "\n")

            #Assign the detected min outlier to the output list
            output_list[[1]] <- lower_outlier_ident
            output_list[[2]] <- lower_outlier_index

            #Name the output list elements
            names(output_list) <- c("left-tail outlier", "element index")

        } else {
            message("No lower outlier have been identified. \n")

        }


    } else {
        message("Identifying the highest outlier on the right tail...")
        
        #An if statement to determine if there is an outlier
        if (max_z_score > crit_t) {
            
            #Assign the outlier identity and index
            upper_outlier_ident <- non_normal_list_element$f_constants[z_scores %in% max_z_score]
            upper_outlier_index <- match(x = upper_outlier_ident, non_normal_list_element$f_constants)

            message("The following maximum value was found to be an outlier: ", upper_outlier_ident, "\n", "at the following row index: ", upper_outlier_index, "\n")
            
            #Assign the detected min outlier to the output list
            output_list[[1]] <- upper_outlier_ident
            output_list[[2]] <- upper_outlier_index

            #Name the output list elements
            names(output_list) <- c("right-tail outlier", "element index")

        } else {
            message("No upper outlier have been identified. \n")

        }


    }


    ##Return the output list
    return(output_list)

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
    f_constants <- normal_list_element$f_constants

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
    f_constants <- normal_list_element$f_constants

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
                
            }


        } else {
            message("No significant outlier found on the left tail.")

        }

    } else {
        
        #Assign the value of Grubbs test on the right tail
        grubbs_upper_res <- outliers::grubbs.test(f_constants)

        #This if statement controls that if there is an outlier it should be identified and removed
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
                
            }

        } else {
            message("No significant outlier found on the right tail.")
            
        }
    }

    #Return the output list
    return(output_list)
    
}




## This a wrapper function which binds the 3 outlier detection and removal functions together and runs them in a  while loop until there are no more
## outliers
## IMPORTANT NOTE: this function seems to work, but is a bit cobbled together as I'm not very familiar with while loops. Therefore please if you have more
## experience with them check and debug this bad boy! :)
outlier_cleaner = function(calculate_f_constants_output_list, is_data_normal_output_list, nonparam_test = "numeric_outlier_test") {
    for (index in seq_along(is_data_normal_output_list)) {
        if (is_data_normal_output_list[[index]]$p.value < 0.05) {
            message("It seems your f-constants show a non-normal distribution. The chosen non-parametric test will be used for outlier detection and removal. \n")

            if (nonparam_test == "numeric_outlier_test") {
                for (e in seq_along(calculate_f_constants_output_list)) {
                    calculate_f_constants_output_list[[e]] <- remove_outlier_f_const_NOTest(calculate_f_constants_output_list[[e]])

                }

                return(calculate_f_constants_output_list)
                
            } else {
                while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- remove_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[i]], left_tail = TRUE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

                while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- remove_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[i]], left_tail = FALSE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

                return(calculate_f_constants_output_list)
            }


        } else {
            message("It seems your f-constants show a normal distribution. A parametric Grubbs test will be used for outlier detection and removal. \n")

            while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- remove_outlier_f_const_Grubbs(calculate_f_constants_output_list[[i]], left_tail = TRUE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

                while_loop_run <- TRUE
                while (while_loop_run == TRUE) {
                    for (i in seq_along(calculate_f_constants_output_list)) {
                        result <- remove_outlier_f_const_Grubbs(calculate_f_constants_output_list[[i]], left_tail = FALSE)
                        calculate_f_constants_output_list[[i]] <- result[[1]]
                        while_loop_run <- result[[2]]  # Update while_loop_run based on outliers_found
                    }
                }

            return(calculate_f_constants_output_list)
        }
    }
}




## This a wrapper function which binds the 3 outlier detection functions together
outlier_detector = function(calculate_f_constants_output_list, is_data_normal_output_list, nonparam_test = "numeric_outlier_test", left_tail = TRUE) {
    ##Declare variables

    #Declare the return list for the numeric outlier test (NOTest)
    return_list_NOTest <- vector(mode = "list", length = 2)
    
    #Declare the return lists for the modified Z-score test
    return_list_mZscore <- vector(mode = "list")

    #Declare the return lists for the Grubbs test
    return_list_Grubbs <- vector(mode = "list")


    ##Run the appropriate detector functions

    #This construct checks if the f-constants from the elements of the calculate f-constants output list are normally distributed
    #NOTE: this is possible because the element order is preserved between the calc.f-const list and the is_data_normal list 
    for (index in seq_along(is_data_normal_output_list)) {
        
        #This if statement checks if the f-factor distribution is normal
        if (is_data_normal_output_list[[index]]$p.value < 0.05) {
            message("It seems your f-constants show a non-normal distribution. The chosen non-parametric test will be used for outlier detection. \n")
        
            #This if statement is responsible for choosing a non-parametric test
            if (nonparam_test == "numeric_outlier_test") {
                for (e in seq_along(calculate_f_constants_output_list)) {
                    return_list_NOTest[[e]] <- detect_outlier_f_const_NOTest(calculate_f_constants_output_list[[e]])

                }

                #Return the detected outlier identities and indexes
                return(return_list_NOTest)
                
            } else {

                #Doing the outlier test on the left tail with the modified Z-score test
                for (i in seq_along(calculate_f_constants_output_list)) {
                    return_list_mZscore[[i]] <- detect_outlier_f_const_mZscore_test(calculate_f_constants_output_list[[i]], left_tail = left_tail)

                }

                #Return the detected outlier identities and indexes
                return(return_list_mZscore)

            }
        
        } else {
            message("It seems your f-constants show a normal distribution. A parametric Grubbs test will be used for outlier detection. \n")

            #Doing the outlier test on the left tail with Grubbs test
            for (i in seq_along(calculate_f_constants_output_list)) {
                return_list_Grubbs[[i]] <- detect_outlier_f_const_Grubbs(calculate_f_constants_output_list[[i]], left_tail = left_tail)
                        
            }

            #Return the detected outlier identities and indexes
            return(return_list_Grubbs)

        }

    }   

}


################################################# Section end #################################################





################################################    MARK: calculate the mean f-constants    #################################################
                                                #       and estimate tumor volumes          #




## Thi function will calculate the mean  f-constants
calc_mean_f = function(calculate_f_constants_output_list) {
    ##Declare the function variables

    #List to modify
    input_list <- vector(mode = "list", length = length(calculate_f_constants_output_list))

    #List to return
    return_list <- vector(mode = "list", length = length(calculate_f_constants_output_list))


    ##Assign function variables
    input_list <- calculate_f_constants_output_list


    ##Calculations
    for (index in seq_along(input_list)) {
        return_list[[index]] <- mean(input_list[[index]]$f_constants)
    }


    ##Return the output list
    return(return_list)

}




## This function will estimate the tumor volumes based on the mean f-constant
estimate_tumor_volume = function(calculate_f_constants_output_list, mean_f_values_list) {
    ##Declare the function variables

    #List to modify
    input_list <- vector(mode = "list", length = length(calculate_f_constants_output_list))

    #List to return
    return_list <- vector(mode = "list", length = length(calculate_f_constants_output_list))


    ##Assign function variables
    input_list <- calculate_f_constants_output_list


    ##Calculations
    
    #This for loop will walk along the input list to access each element
    #NOTE: the length of the input_df and the mean_f_values_list is the same, so indexes can be used for both
    for (index in seq_along(input_list)) {
        #Declare a temporary dataframe
        temp_df <- data.frame(matrix(nrow = nrow(input_list[[index]]), ncol = ncol(input_list[[index]])))

        #Assign temp_df
        temp_df <- input_list[[index]]

        #This for loop will progress through the rows of each individual list element (dataframe) to estimate the tumor volumes using the 
        #appropriate formula
        for (ind in seq_len(nrow(input_list[[index]]))) {
            
            #Estimate and assign the tumor volumes to each individual sample in the temp_df
            temp_df$Estim_tumor_vol[ind] <- (pi / 6) * mean_f_values_list[[index]] * (temp_df$L[ind] * temp_df$W[ind])^(3 / 2)
        }

        #Assign the temp_df to the return list for each input list element
        return_list[[index]] <- temp_df

    }


    ## Return the modified dataframes
    return(return_list)

}


################################################# Section end #################################################





################################################    MARK: correct the estimated volumes    #################################################
                                                #       and test the goodness of fit       #




## Correct the estimated tumor volumes
## NOTE: this is a complex function. The idea is that it tests what percentage of tumor volume deviation comes from the 0.5, 1, 1.5, 2 
## standard deviation distance of sample f-constants from the calculated mean f-constant. Then based on which SD bracket
##f-constants fall into it corrects the estimated tumor volume by the appropriate percentage depending on if the value
## is under or over estimated (this depends on which direction the f-constant is deviating from the mean)
correct_tumor_vol = function(estimated_tumor_volume_list, mean_f_values_list) {
    ##Declare static function variables

    #Define the input list
    input_list <- vector(mode = "list", length = length(estimated_tumor_volume_list))

    #Define the return/output list
    return_list <- vector(mode = "list", length = length(estimated_tumor_volume_list))


    ##Assign the appropriate variables

    #Assign the input list
    input_list <- estimated_tumor_volume_list


    ##Calculate the appropriate variables

    for (index in seq_along(input_list)) {
        ##Declare dynamic function variables

        #Define a vector for the standard deviation fractions
        std_frac <- vector(mode = "numeric", length = 5)

        #Define a vector containing the standard deviations of 0.5-2 SD below mean
        below_mean_std_vals <- vector(mode = "numeric", length = 5)
        names(below_mean_std_vals) <- as.character(seq(from = 0, to = 2, by = 0.5))

        #Define a vector containing the standard deviations of 0.5-2 SD above mean
        above_mean_std_vals <- vector(mode = "numeric", length = 5)
        names(above_mean_std_vals) <- as.character(seq(from = 0, to = 2, by = 0.5))

        #Define a vector containing the tumor volume deviation (in percentage) of 0.5-2 SD below mean
        below_mean_volumes <- vector(mode = "numeric", length = 5)
        names(below_mean_volumes) <- as.character(seq(from = 0, to = 2, by = 0.5))

        #Define a vector containing the tumor volume deviation (in percentage) of 0.5-2 SD above mean
        above_mean_volumes <- vector(mode = "numeric", length = 5)
        names(above_mean_volumes) <- as.character(seq(from = 0, to = 2, by = 0.5))

        #Define the temp_df
        temp_df <- data.frame(matrix(nrow = nrow(input_list[[index]]), ncol = ncol(input_list[[index]])))

        #Define a dynamic filter_df
        filter_list <- vector(mode = "list", length = 4)

        #Define the standard deviation variable
        temp_std <- vector(mode = "numeric", length = 1)
        

        ##Assign the dynamic function variables

        #Assign the temp_df 
        temp_df <- input_list[[index]]

        #Assign the standard deviation variable
        temp_std <- sd(temp_df$f_constants)

        #Assign the standard deviation fractions
        std_frac <- seq(from = 0, to = 2, by = 0.5)

        #Assign the standard deviations below mean (we start with 0 - the mean - and go down by half a sd compared to the mean)
        below_mean_std_vals <- mean_f_values_list[[index]] - (temp_std * std_frac)

        #Assign the standard deviations above mean (we start with 0 - the mean - and go up by half a sd compared to the mean)
        above_mean_std_vals <- mean_f_values_list[[index]] + (temp_std * std_frac)

        #Assign the tumor volume deviation below mean (in this case the tumor volumes are
        #over-estimated)
        for (ind in seq_along(below_mean_std_vals)) {
            filter_list[[ind]] <- dplyr::filter(temp_df, f_constants < below_mean_std_vals[ind] & f_constants > below_mean_std_vals[ind + 1])
            
            below_mean_volumes[ind] <- sum(filter_list[[ind]]$Estim_tumor_vol) / sum(filter_list[[ind]]$Tumor_volume_.mm3.)
            
        }

        #Assign the tumor volume deviation above mean (in this case the tumor volumes are
        #under-estimated)
        for (ind in seq_along(above_mean_std_vals)) {
            filter_list[[ind]] <- dplyr::filter(temp_df, f_constants > above_mean_std_vals[ind] & f_constants < above_mean_std_vals[ind + 1])
            
            above_mean_volumes[ind] <- sum(filter_list[[ind]]$Estim_tumor_vol) / sum(filter_list[[ind]]$Tumor_volume_.mm3.)
            
        }

        #Perform the value correction for values which f-constants are below the mean f-constant (in this case the tumor volumes are
        #over-estimated). The outer for loop traverses the rows of the temp_df
        for (i in seq_len(nrow(temp_df))) {
            
            #The inner for-loop traverses the standard deviation values of f-constants below the mean f-constant (the SD here is calculated for the 
            #f-constants and brackets of 0(mean)-0.5-1-1.5-2 SD are used)
            for (e in seq_along(below_mean_std_vals)) {
                
                #An if statement to control which correction factor (his is volume parentage of the estimated volume compared to the measured volume)
                #should be used, based on which f-constant bracket the sample falls into
                if (temp_df$f_constants[i] < below_mean_std_vals[e] && temp_df$f_constants[i] > below_mean_std_vals[e + 1]) {
                    
                    #Perform the correction by dividing the estimated tumor volume with the correction factor (the percentage of over or under estimation
                    #as a function of the f-constant deviation) and assign the corrected values into a new column
                    temp_df$Corrected_tumor_vol[i] <- temp_df$Estim_tumor_vol[i] / below_mean_volumes[e]

                    #assign the correction factors into a new column
                    temp_df$Correction_factor[i] <- below_mean_volumes[e]
                }

            }
            
        }

        #Perform the value correction for values which f-constants are above the mean f-constant (in this case the tumor volumes are
        #under-estimated). The outer for loop traverses the rows of the temp_df
        for (i in seq_len(nrow(temp_df))) {
            
            #The inner for-loop traverses the standard deviation values of f-constants above the mean f-constant (the SD here is calculated for the 
            #f-constants and brackets of 0(mean)-0.5-1-1.5-2 SD are used)
            for (e in seq_along(above_mean_std_vals)) {
                
                #An if statement to control which correction factor (his is volume parentage of the estimated volume compared to the measured volume)
                #should be used, based on which f-constant bracket the sample falls into
                if (temp_df$f_constants[i] > above_mean_std_vals[e] && temp_df$f_constants[i] < above_mean_std_vals[e + 1]) {
                    
                    #Perform the correction by dividing the estimated tumor volume with the correction factor (the percentage of over or under estimation
                    #as a function of the f-constant deviation) and assign the corrected values into a new column
                    temp_df$Corrected_tumor_vol[i] <- temp_df$Estim_tumor_vol[i] / above_mean_volumes[e]

                    #assign the correction factors into a new column
                    temp_df$Correction_factor[i] <- above_mean_volumes[e]
                }

            }
            
        }

        
        ##Add the temporary, modified dataframes to the return list
        return_list[[index]] <- temp_df
        

    }

    
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
