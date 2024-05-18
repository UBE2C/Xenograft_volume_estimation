#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#     COMIX data processer module       #
#                                       #
#########################################






################################################# MARK: Package and path #################################################
                                                #       management       #





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




## Define the paths for the script and I/O files


# Define the path to the directory this script is in
script_dir_path <- this.path::here()
setwd(script_dir_path)

# Define the preferred input directory path
pref_input_path <- paste0(script_dir_path, "/", "Input_files")

# Define the preferred final output directory path
pref_output_path <- paste0(script_dir_path, "/", "Intermediate_IO")

# Define the preferred path to the intermediate input/output file directory for later use
final_output_path <- paste0(script_dir_path, "/", "Output_files")



# Define the input directory path
def_input_path <- script_dir_path

# Define the final output directory path
def_output_path <- script_dir_path


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
        # Define the input directory path
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



## Define the default output file name
def_output_name <- "output_table"


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




## Add command line arguments


# Create the arguments list
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
    By default the option is set to FALSE.")

)


# Create a program description
prog_descr <- c(paste0("The COMIX-data_processer is specifically made to process the caliper tumor measurements done for the COMIX experiments.\n",
                "This program can be run form the command line as a standalone program or as a subprogram for the full data processing and tumor volume calculation script. \n",
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
                                                #   and processing   #




## Data reading function
read_data = function(csv_file_path, verbose) {
    
    #List files in the file path
    csv_files <- list.files(csv_file_path, full.names = TRUE)
    
    #file_list <- vector(mode = "character")
    input_list <- list()

    #A for loop to handle multiple input files
    #NOTE: the loop will assign each .csv file into a new dataframe in the global env
    for (index in seq_along(csv_files)) {
        #An if statement to control non .csv inputs
        if (grepl(".csv$", csv_files[index]) == TRUE) {
            tmp_file <- read.csv(file = csv_files[index], header = TRUE, sep = ";")
            input_list[[index]] <- tmp_file
            
        } else {
            if (verbose == TRUE) {
                warning(paste0("The input file ", csv_files[index], " is not a .csv.", "\n",
                "Jumping over to the next file."))
                next
            
            } else {
                next
            
            } 
            
            
        }
    
        
    }

    #An sapply subset to remove NULL elements form the list (these were created from non .csv entries)
    output_csv_list <- input_list[!sapply(input_list, is.null)]
    
    #An if statement controlling how much information the function should return if verbose is TRUE
    if (verbose == TRUE) {
        message("\nPrinting the loaded dataframes ready for processing:")
        
        for (e in seq_along(output_csv_list)) {
            print(output_csv_list[[e]])
        }
        
        
    }

    #Return the output list
    return(output_csv_list)

}




##  This function convert the L x W notation into the actual multiplication value
process_data = function(data_list, verbose) {
    #Create a new dataframe with the sample IDs and L-W entries
    processed_data <- data.frame()
    processed_data_lst <- list()

    #This loop will load the dataframes to process one by one
    for (index in seq_along(data_list)) {
        tmp_dataframe <- data_list[[index]]
        tmp_list <- list()

        #Determining the range for the second for loop based on the names of the measurement columns
        measurement_cols <- grep(pattern = "^X", x = colnames(tmp_dataframe))

        #This for loop will go though the appropriate columns and will split the entries and adds them to
        #an empty dataframe
        for (element in measurement_cols) {
                
                #This will does the L x W multiplication and assigns the resulting vectors into the temporary list
                tmp_list[[element]] <-  as.numeric(stringr::str_split_i(string = tmp_dataframe[, element],
                pattern = "x", i = 1)) * as.numeric(stringr::str_split_i(string = tmp_dataframe[, element],
                pattern = "x", i = 2))

                #convert the list tmp_list elements into a data frame and store it in the processed_data df
                processed_data <- as.data.frame(do.call(cbind, tmp_list))
            
        }

        #Add the processed_data dfs to the processed list and set the names of each column to the original col name
        processed_data_lst[[index]] <- processed_data
        colnames(processed_data_lst[[index]]) <- colnames(tmp_dataframe)[measurement_cols]

        #Add the rownames as a new column so that they will be saved when the data is written into new .csv files
        processed_data_lst[[index]] <- dplyr::mutate(.data = processed_data_lst[[index]],
                                        Mouse_ID = data_list[[index]]$Mouse.ID,
                                        Treatment_group_ID = paste0("G", index),
                                        Treatment = data_list[[index]]$Treatment,
                                        .before = 1)

        #Modify the column names to reflect that the L and W values are multiplied together
        colnames(processed_data_lst[[index]]) <- stringr::str_replace_all(colnames(processed_data_lst[[index]]), "^X", "LxW_")

    }
    
    #An if statement controlling how much information the function should return if verbose is TRUE
    if (verbose == TRUE) {
        message("\nPrinting the processed dataframes:")
        for (index in seq_along(processed_data_lst)){
            print(processed_data_lst[[index]])
        }
    }

    #Return the list containing the processed DFs
    return(processed_data_lst)

}


################################################# Section end #################################################






################################################# MARK: Main function #################################################
                                                #                     #




## The main function which will return the clean, modified data frame/data frames
main = function(input_path = arguments$input_path, output_path = arguments$output_path, 
                output_name = arguments$output_name, verb = arguments$verbose) {

    #Call the package management functions with an option to verbose or not
    if (verb == TRUE) {
        package_controller(CRAN_packages)
        package_loader(CRAN_packages)        
    
    } else {
        suppressMessages(package_controller(CRAN_packages))
        suppressPackageStartupMessages(package_loader(CRAN_packages))
    }

    #Call the read_data function to load the data from the default or given input path
    input_csv_list <- read_data(input_path, verbose = verb)
    print(input_csv_list)

    #Call the processing function with the output
    processed_csv_list <- process_data(input_csv_list, verbose = verb)
    print(processed_csv_list)
    
    #Save the processed files based on the output path
    if (dir.exists(output_path) == TRUE) {
        for (index in seq_len(length(processed_csv_list))) {
            readr::write_csv(x = processed_csv_list[[index]], file = paste0(output_path, "/", "processed", "_", output_name, "_", index, ".csv"))
        }

    } else {
        #An if statement controlling how much information the function should return if verbose is TRUE
        if (verb == TRUE) {
            message("\nThe given directory on the output_path does not exists. Creating one...")
        }

        #Creating the given directory as a save location
        dir.create(output_path)
        
        #Saving the processed .csv files to the output directory
        for (index in seq_len(length(processed_csv_list))) {
            readr::write_csv(x = processed_csv_list[[index]], file = paste0(output_path, "/", output_name, "_", index, ".csv"))
        }

    }

    #Call the session info if verbose is TRUE
    if (verb == TRUE) {
        message("\nPrinting the R session info:")
        print(sessionInfo())
    }

    message(paste0("\nThe data processing was successful! \n",
                    "Your processed data was saved to:\n",
                    output_path, "\n",
                    "with the name: ", "processed", "_", output_name, "\n"))
    
}


################################################# Section end #################################################






################################################# MARK: call the Main function #################################################
                                                #                              #




## Calling the main function
main()


################################################# Section end #################################################
