#!user/bin/env -S Rscript --vanilla


#########################################
#                                       #
#   Tumor measurement data processing   #
#               COMIX data              #
#                                       #
#########################################



## Check for the required packages and install them if missing


# The list of required packages
packages <- c("tidyverse", "optparse")


# A function to check and install packages
package_controller = function(packages){
    
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
                message("The requested package is already installed.")
            }
        }
    } else {
        message("No package was requested.")
    }
    
}


# A function to load required packages
package_loader = function(packages) {
    lapply(packages, library, character.only = TRUE)
}




## Define the path to the directory this script is in
script_path <- rstudioapi::getActiveDocumentContext()$path
script_dir_path <- stringr::str_remove(path, "[A-z0-9]+\\W[A-z0-9]+\\W[A-z0-9]+.r$")
setwd(script_dir_path)




## Data reading function
read_data = function(csv_file_path = script_dir_path) {
    
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
            warning(paste0("The input file ", csv_files[index], " is not a .csv.", "\n",
            "Jumping over to the next file."))
            next
        }
    
        
    }

    #An sapply subset to remove NULL elements form the list (these were created from non .csv entries)
    output_csv_list <- input_list[!sapply(input_list, is.null)]
    
    #Return the output list
    return(output_csv_list)
}




##  This function splits the Length and Width measurements in the loaded tumor measurement files
process_data = function(data_list) {
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

    #Return the list containing the processed DFs
    return(processed_data_lst)
}




## The main function which will return the clean, modified data frame/data frames
main = function() {
    args <- commandArgs(trailingOnly = TRUE)

    # Call the package management functions
    package_controller()
    package_loader()
    



    
}








