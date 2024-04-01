#!user/bin/Rscript


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




## Data reading function
read_data = function(csv_file_path) {
    
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
            #assign(x = csv_files[index], value = tmp_file, envir = .GlobalEnv)

            #Creating a file list without NA values as an input for the process function
            #file_list[index] <- csv_files[index]
            #file_list <- file_list[!is.na(file_list)]
            #assign(x = "file_list", value = file_list, envir = .GlobalEnv)
        } else {
            warning(paste0("The input file ", csv_files[index], " is not a .csv.", "\n",
            "Jumping over to the next file."))
            next
        }
    
        
    }

    #An sapply subselect to remove NULL elements form the list (these were created from non .csv entries)
    input_list <- input_list[!sapply(input_list, is.null)]
    
    #Return the output list
    return(input_list)
}




##  This function splits the Length and Width measurements in the loaded tumor measurement files
process_data = function(file_list = file_list) {
    #Create a new dataframe with the sample IDs and L-W entries
    processed_data <- data.frame(matrix(nrow = 14, ncol = 22))
    test_lst <- list()

    #This loop will load the dataframes to process one by one
    for (index in seq_along(file_list)) {
        tmp_dataframe <- get(file_list[index])

        #Determining the range for the second for loop based on the names of the measurement columns
        measurement_cols <- grep(pattern = "^X", x = colnames(tmp_dataframe))

        #This for loop and if statement will go though the appropriate columns and will split the entries
        for (i in seq_along(measurement_cols)) {
            for (e in seq_along(nrow(processed_data))){
                processed_data[e, ] <- stringr::str_split_i(string = tmp_dataframe[, measurement_cols[i]],
                pattern = "x", i = 1)

                processed_data[e + 1, ] <- stringr::str_split_i(string = tmp_dataframe[, measurement_cols[i]],
                pattern = "x", i = 2)
            }
            #processed_data[index][[i]] <- rbind(stringr::str_split(string = tmp_dataframe[, measurement_cols[i]],
            #pattern = "x", simplify = TRUE))
            
            #print(cbind(stringr::str_split(string = tmp_dataframe[, measurement_cols[i]],
            #pattern = "x", simplify = TRUE)))
            
        }
    
    test_lst[[index]] <- processed_data
    print(test_lst)
    #Return the processed data
    #return(processed_data)

    }

    
}



##  This function splits the Length and Width measurements in the loaded tumor measurement files
process_data = function(data_list) {
    #Create a new dataframe with the sample IDs and L-W entries
    processed_data <- data.frame()
    #test_lst <- list()

    #This loop will load the dataframes to process one by one
    for (index in seq_along(data_list)) {
        tmp_dataframe <- get(data_list[index])

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

    #This snippet will split up the resulting unified dataframe to the two input files
    processed_group_1 <- processed_data[1:14, ]
    processed_group_2 <- processed_data[15:28, ]

    #Prepare the proper column and row names for Group1
    split_measurements <- c("_L", "_W")
    g1_sample_dates <- colnames(get(data_list[1])[7:13])
    
    g1_rownames <- paste0(rep(g1_sample_dates, each = 2), rep(split_measurements, times = 7))
    g1_colnames <- get(data_list[1])[, 3]

    #Prepare the proper column and row names for Group2
    g2_sample_dates <- colnames(get(data_list[2])[7:13])
    
    g2_rownames <- paste0(rep(g2_sample_dates, each = 2), rep(split_measurements, times = 7))
    g2_colnames <- get(data_list[2])[, 3]

    #Assign the new column and row names to the new dataframes
    rownames(processed_group_1) <- g1_rownames
    colnames(processed_group_1) <- g1_colnames

    rownames(processed_group_2) <- g2_rownames
    colnames(processed_group_2) <- g2_colnames

#################################################################################################################################
#continue from here
    print(processed_group_1)
    print(processed_group_2)
}




## The main function which will return the clean, modified data frame/data frames
main = function() {
    args <- commandArgs(trailingOnly = TRUE)

    # Call the package management functions
    package_controller()
    package_loader()
    


    
}








