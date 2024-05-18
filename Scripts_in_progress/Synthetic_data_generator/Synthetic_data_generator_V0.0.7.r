#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#       Synthetic data generation       #
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






################################################# MARK: Subfunction - Dates #################################################
                                                #                           #




## This subfunction creates a series of requested dates as  colnames for the synthetic tumor volume data
## NOTE: the function returns a vector of dates (type numeric)
create_dates = function(no_of_days, start_date = "2024-01-01", by_days) {
    dates <- vector(mode = "numeric", length = no_of_days)

    dates <- seq(from = as.Date(start_date), length.out = no_of_days, by = by_days)

    return(dates)
}


################################################# Section end #################################################






################################################# MARK: Subfunction - Data_generator #################################################
                                                #                                    #




## This is the main subfunction responsible for creating synthetic tumor volume measurements. The function is supposed to model the growth of
## subcutaneously injected xenograft tumors under various condition (control vs treatment)
## NOTE: the function returns a list with two dataframes
tumor_data_generator = function(number_of_measurements, number_of_samples, request_dates = FALSE, dates, seed = 1234, initial_element_mean = 1.5,
    initial_element_sd = 0.3, growth_variation_min = 0.2, growth_variation_max = 0.5, treatment_administered = FALSE, treatment_effect_start = 3, 
    treatment_variation_min = -0.1, treatment_variation_max = 0.1, continuous_treatment = TRUE) {
    
    ##Declare function variables

    #A vector containing the factor for tumor growth variation
    growth_variation <- vector(mode = "numeric", length = 1)

    #A vector containing the factor for the treatment effect
    treatment_effect <- vector(mode = "numeric", length = 1)

    #A vector containing the iteration the treatment effect starts
    treatment_effect_iteration <- vector(mode = "numeric", length = 1)

    #A vector containing the number of measurements (reflects caliper measurements)
    nMeasurements <- vector(mode = "numeric", length = 1)
    
    #A vector containing the names of measurement dimensions (L, W, H)
    dimensions <- vector(mode = "character", length = 3)

    #A vector containing the "caliper measurements" 
    measurement <- vector(mode = "numeric", length = number_of_measurements)

    #A vector containing the number of samples measured
    nSamples <- vector(mode = "numeric", length = 1)
    
    #A list representing a single sample, containing the corresponding L, W and H measurements
    sample <- list()

    #A list containing the calculated volumes for each sample
    sample_volumes_list <- vector(mode = "list", length = number_of_measurements)

    #A list containing the LxW variables for the f-factor estimation
    LxW_values_list <- vector(mode = "list", length = number_of_measurements)

    #Vectors containing the names of the samples (rownames) and the measurement dates (colnames) for the output dataframe
    sample_names <- vector(mode = "character", length = number_of_samples)
    measurement_dates <- vector(mode = "character", length = length(number_of_measurements))

    #The output dataframes containing the generated synthetic tumor volume or LxW measurements
    output_volume_dataframe <- data.frame(matrix(nrow = number_of_samples, ncol = length(number_of_measurements)))
    output_LxW_dataframe <- data.frame(matrix(nrow = number_of_samples, ncol = length(number_of_measurements)))

    #A list containing the final output dataframes
    final_output_list <- vector(mode = "list", length = 2)
    
    
    ##Initialize function variables

    #Initialize the treatment effect iteration variable
    treatment_effect_iteration <- treatment_effect_start
    
    #Initialize the number of measurements variable
    nMeasurements <- number_of_measurements

    #initialize the names of the measurement dimensions variable
    dimensions <- c("L", "W", "H")

    #Initialize the number of samples variable
    nSamples <- number_of_samples
    
    #Initialize the sample_names variable with a for-loop
    for (number in seq_len(nSamples)) {
        sample_names[number] <- paste0("S", number)
    }

    #Initialize the measurement dates variable if requested
    if (request_dates == TRUE) {
        measurement_dates <- dates
    }
    
    
    
    ##Set random number generation seed
    set.seed(seed = seed)

    
    ##Generate sample data
    
    #The main for loop which controls the number of samples (iterations) needed
    for (element in seq_len(nSamples)) {
        #Declare the dynamic volume and LxW variables (this needs to be reset for every sample, hence it is declared inside the main loop)
        volume <- vector(mode = "numeric", length = nMeasurements)
        LxW_value <- vector(mode = "numeric", length = nMeasurements)

        #This inner loop is responsible for creating a single sample which has the corresponding number of L, W and H measurements (dimensions)
        for (dim in dimensions) {
            
            #Declare the dynamic initial_element variable (this needs to be reset for every dimension, hence it is declared inside the loop)
            initial_element <- vector(mode = "numeric", length = 1)

            #This while loop ensures tha the initial element is more than or equal to 0.5 (to model some level of tumor growth to the first measurement
            #and to prevent a negative tumor dimension)
            while (initial_element <= 0.5) {
            initial_element <- rnorm(n = 1, mean = initial_element_mean, sd = initial_element_sd)
            }

            #This is the innermost for loop responsible for generating the random measurements for a single dimension (either L, W or H)
            for (index in seq_len(nMeasurements - 1)){
                
                #Declare the dynamic growth_variation variable (this needs to be reset for every measurement of a dimension, hence it is declared inside the loop)
                growth_variation <- runif(n = 1, min = growth_variation_min, max = growth_variation_max)
                
                #Create the first measurement of the dimension vector
                measurement[1] <- initial_element

                #Add the rest of the measurements which mimic the growth of the tumor by a random % range (the % range is given as the growth_variation min/max values)
                measurement[index + 1] <- measurement[index] + measurement[index] * growth_variation
            
            }

            #An if statement controlling if a treatment should be taken into consideration (the treatment should slow down the tumor growth)
            if (treatment_administered == TRUE) {
                
                #An if statement controlling if a continuous treatment administration should be taken into consideration or just a single dose
                if (continuous_treatment == TRUE) {

                    #This is the innermost for loop responsible for generating the random measurements for a single dimension (either L, W or H)
                    #weighted by the treatment effect
                    for (index in seq(from = treatment_effect_iteration, to = nMeasurements, by = 1)) {
                    
                        #Declare the dynamic growth_variation variable (this needs to be reset for every measurement of a dimension, hence it is declared inside the loop)
                        treatment_effect <- runif(n = 1, min = treatment_variation_min, max = treatment_variation_max)

                        #Set the element where the treatment was administered to the new initial element
                        measurement[treatment_effect_iteration] <- measurement[treatment_effect_iteration] + measurement[treatment_effect_iteration] * treatment_effect

                        #Create the rest of the measurements of the given dimension to mimic the treatment effect by a random % range
                        measurement[index + 1] <- measurement[index] + abs(measurement[index]) * treatment_effect

                    }

                } else {
                    
                    #This is the innermost for loop responsible for modifying the previously generated measurements for a single dimension (either L, W or H)
                    for (index in seq(from = treatment_effect_iteration, to = nMeasurements, by = 1)) {
                        
                        #Declare the dynamic growth_variation variable (this needs to be reset for every measurement of a dimension, hence it is declared inside the loop)
                        treatment_effect <- runif(n = 1, min = treatment_variation_min, max = treatment_variation_max)

                        #Modify the rest of the measurements of the given dimension to mimic the treatment effect by a random % range
                        measurement[index] <- measurement[index] + measurement[index] * treatment_effect

                    }

                }

                #A quality control loop responsible to set possible negative growth measurements to 0
                for (e in seq_along(measurement)) {
                    if (measurement[e] < 0) {
                        measurement[e] <- 0
                    }
                }

            }

            #Add the measurements to the L, W, H dimensions of a sample
            sample[[dim]] <- measurement
            
        }

        #Calculate the volume of every sample in a measurement and add them to the volume variable (each element is a tumor volume for a given date/measurement)
        for (ind in seq_len(nMeasurements)) {
            volume[ind] <- (4 / 3) * pi * (sample$L[ind] / 2) * (sample$W[ind] / 2) * (sample$H[ind] / 2)
            LxW_value[ind] <- (sample$L[ind]) * (sample$W[ind])

        }
        
        #Add the calculated tumor volume measurements to a list of samples
        sample_volumes_list[[element]] <- volume

        #Add the calculated LxW values to a list of samples
        LxW_values_list[[element]] <- LxW_value

    }

    #Convert the output volume list into a volume dataframe
    output_volume_dataframe <- as.data.frame(do.call(rbind, sample_volumes_list))

    #Convert the output LxW list into an LxW dataframe
    output_LxW_dataframe <- as.data.frame(do.call(rbind, LxW_values_list))

    #An if statement controlling if the cols should be named with dates
    if (request_dates == TRUE) {
        #Name the rows and columns of the output dataframes
        colnames(output_volume_dataframe) <- measurement_dates
        rownames(output_volume_dataframe) <- sample_names
        colnames(output_LxW_dataframe) <- measurement_dates
        rownames(output_LxW_dataframe) <- sample_names

    } else {
        #Name the rows of the output dataframes
        rownames(output_volume_dataframe) <- sample_names
        rownames(output_LxW_dataframe) <- sample_names

    }
    

    #Add the output dataframe to the final output list
    final_output_list[[1]] <- output_volume_dataframe
    final_output_list[[2]] <- output_LxW_dataframe

    #Name the final output list elements
    names(final_output_list) <- c("Tumor_volumes", "Tumor_LxW_values")
    

    ##Return the output dataframe
    return(final_output_list)

}


################################################# Section end #################################################






################################################# MARK: Subfunctions - Data_formatters #################################################
                                                #                                    #





## This function formats the Tumor_volume and Tumor_LxW_values dataframes to mimic the uCT- and caliper-measurements dataframes of the f-factor estimating module
format_tumor_data = function(data, treatment_group_id = "Ctrl", treatment = "control") {
    ##Declare function variables

    #Declare the output list
    output_list <- vector(mode = "list", length = 2)


    for (i in seq_along(data)) {
        ##Declare dynamic function variables
        

        #Declare and initialize the input dataframe containing the volume measurements for transformation
        input_data <- data[[i]]

        #Declare a vector for the modified column names of the transformed df
        new_column_names <- vector(mode = "character", length = ncol(input_data))

        #Declare the final output dataframe 
        output_data <- data.frame(matrix(nrow = nrow(input_data), ncol = ncol(input_data) + 3))


        ##Initialize function variables

        #Initialize the sample IDs to be added to the transformed df
        sample_id <- rownames(input_data)
        
        if (names(data[i]) == "Tumor_volumes") {
            #Initialize the modified column names of the transformed df
            new_column_names <- paste0("uCT_volume", "_", colnames(input_data))
        } else {
            #Initialize the modified column names of the transformed df
            new_column_names <- paste0("LxW", "_", colnames(input_data))

        }
        


        ##Transform the input dataframe

        #Add the new column names to the input dataframe
        colnames(input_data) <- new_column_names

        #Transform the input dataframe by adding additional columns to it
        output_data <- dplyr::mutate(.data = input_data, Treatment_group_id = treatment_group_id, Treatment = treatment, Mouse_ID = sample_id, .before = 1)

        #Add the transformed dataframe to the output list
        output_list[[i]] <- output_data

    }
    

    ##Return the transformed dataframes
    return(output_list)

}


################################################# Section end #################################################






################################################# MARK: Main function #################################################
                                                #                                    #





##The main function integrating all of the subfunctions
main = function(number_of_measurements, number_of_samples, request_dates = FALSE, dates, seed = 1234, initial_element_mean = 1.5,
    initial_element_sd = 0.3, growth_variation_min = 0.2, growth_variation_max = 0.5, treatment_administered = FALSE, treatment_effect_start = 3, 
    treatment_variation_min = -0.1, treatment_variation_max = 0.1, continuous_treatment = TRUE,
    start_date = "2024-01-01", by_days,
    treatment_group_id = "Ctrl", treatment = "control") {


}

