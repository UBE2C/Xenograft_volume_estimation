#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#     COMIX data processer module       #
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




## Define the default output path and file name
def_output_path <- paste0(script_dir_path, "/", "Intermediate_IO")
def_output_name <- "output_table"



## Add command line arguments


# Create the arguments list
options_list <- list(
    optparse::make_option(opt_str = c("-p", "--input_path"), action = "store", type = "character", default = script_dir_path,
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
        print(output_csv_list[[1]])
        print(output_csv_list[[2]])
        
    }

    #Return the output list
    return(output_csv_list)

}




##  This function splits the Length and Width measurements in the loaded tumor measurement files
process_data = function(data_list, verbose) {
    #Create a new dataframe with the sample IDs and L-W entries
    processed_data <- data.frame()
    processed_data_lst <- list()

    #This loop will load the dataframes to process one by one
    for (index in seq_along(data_list)) {
        tmp_dataframe <- data_list[[index]]

        #Determining the range for the second for loop based on the names of the measurement columns
        measurement_cols <- grep(pattern = "^X", x = colnames(tmp_dataframe))

        #This for loop will go though the appropriate columns and will split the entries and adds them to
        #an empty dataframe
        for (element in measurement_cols) {
                
                #This will split the first numbers representing the Length measurements
                processed_data <- rbind(processed_data, stringr::str_split_i(string = tmp_dataframe[, element],
                pattern = "x", i = 1))

                #This will split the second numbers representing the Width measurements
                processed_data <- rbind(processed_data, stringr::str_split_i(string = tmp_dataframe[, element],
                pattern = "x", i = 2))
            
        }

    }

    #NOTE: This part makes the function useless as a generic data processing function!
    #for a generic use, the following parts need a re-write

    #This snippet will split up the resulting unified dataframe to the two input files
    processed_group_1 <- processed_data[1:14, ]
    processed_group_2 <- processed_data[15:28, ]

    #Prepare the proper column and row names for Group1
    split_measurements <- c("_L", "_W")

    g1_sample_dates <- colnames(data_list[[1]][7:13])
    
    g1_rownames <- paste0(rep(g1_sample_dates, each = 2), rep(split_measurements, times = 7))
    g1_colnames <- data_list[[1]][, 3]

    #Prepare the proper column and row names for Group2
    g2_sample_dates <- colnames(data_list[[2]][7:13])
    
    g2_rownames <- paste0(rep(g2_sample_dates, each = 2), rep(split_measurements, times = 7))
    g2_colnames <- data_list[[2]][, 3]

    #Assign the new column and row names to the new dataframes
    rownames(processed_group_1) <- g1_rownames
    colnames(processed_group_1) <- g1_colnames

    rownames(processed_group_2) <- g2_rownames
    colnames(processed_group_2) <- g2_colnames

    #Add the two processed DFs to the clean processed data list
    processed_data_lst[[1]] <- processed_group_1
    processed_data_lst[[2]] <- processed_group_2

    #Add the rownames as a new column so that they will be saved when the data is written into new .csv files
    processed_data_lst[[1]] <- dplyr::mutate(.data = processed_data_lst[[1]],
                                        Dates = rownames(processed_data_lst[[1]]),
                                        .before = colnames(processed_data_lst[[1]])[1])
    processed_data_lst[[2]] <- dplyr::mutate(.data = processed_data_lst[[2]],
                                        Dates = rownames(processed_data_lst[[2]]),
                                        .before = colnames(processed_data_lst[[2]])[1])

    #An if statement controlling how much information the function should return if verbose is TRUE
    if (verbose == TRUE) {
        message("\nPrinting the processed dataframes:")
        print(processed_data_lst[[1]])
        print(processed_data_lst[[2]])
    }

    #Return the list containing the processed DFs
    return(processed_data_lst)

}




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

    #Call the processing function with the output
    processed_csv_list <- process_data(input_csv_list, verbose = verb)
    
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
                    "with the name: ", output_name, "\n"))
    
}




## Calling the main function
main()







