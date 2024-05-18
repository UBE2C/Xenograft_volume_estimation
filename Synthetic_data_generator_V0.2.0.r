#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#       Synthetic data generation       #
#         by Gabor Bakos (UBE2C)        #
#                                       #
#########################################






################################################# MARK: Package and path #################################################
                                                #       management       #




## Check for the required packages and install them if missing


# The list of required packages
CRAN_packages <- c("tidyverse", "optparse", "this.path", "ggalt")


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

def_output_path <- script_dir_path


################################################# Section end #################################################






################################################# MARK: CLI arguments #################################################
                                                #                     #




## Add command line arguments
options_list <- list(
    optparse::make_option(opt_str = c("-m", "--number_of_measurements"), action = "store", type = "integer", default = 6,
    help = "This argument takes an integer which defines the number measurements the user wants to model.
    By default the script sets the number of measurements to 6."),

    optparse::make_option(opt_str = c("-s", "--number_of_samples"), action = "store", type = "integer", default = 10,
    help = "This argument takes an integer which defines the number samples the user wants to model.
    By default the script sets the number of samples to 10."),

    optparse::make_option(opt_str = c("--growth_model"), action = "store", type = "character", default = "exponential",
    help = "This argument controls if the model should be based on exponential or logistical growth rates. To choose between the growth models type either 'exponential' or 'logistical'.
    By default the growth model is set to 'exponential'."),

     optparse::make_option(opt_str = c("--carrying_capacity"), action = "store", type = "integer", default = 15,
    help = "This argument takes an integer which defines the carrying capacity (a constant which limits tumor volume growth) for the logistical growth model.
    By default the initial carrying capacity is set to 15."),

    optparse::make_option(opt_str = c("--base_intrinsic_growth_rate"), action = "store", type = "double", default = 1,
    help = "This argument takes a floating point number (double) which defines base intrinsic growth-rate for the logistical growth model. 
    This value then will be modified with a random value to represent sample to sample variation in the tumor growth.
    By default the base_intrinsic_growth_rate is set to 0.1."),

    optparse::make_option(opt_str = c("--request_dates"), action = "store_true", default = FALSE, dest = "dates",
    help = "This argument controls if the generated synthetic dataframe should have dates as column names.
    By default the option is set to FALSE."),

    optparse::make_option(opt_str = c("--start_date"), action = "store", type = "character", default = "2024-01-01",
    help = "This argument controls if dates were requested, when the dates should start. The format of the dates should be Y-M-D.
    By default the starting date is '2024-01-01'."),

    optparse::make_option(opt_str = c("--by_days"), action = "store", type = "integer", default = 4,
    help = "This argument controls if dates were requested, how many days should pass between two measurement dates.
    By default the interval is set to 4 days."),

    optparse::make_option(opt_str = c("--seed"), action = "store", type = "integer", default = 1234,
    help = "This argument takes an integer which defines the seed for random number generation. The use of a set seed is important for result consistency
    By default the script sets the seed number to 1234."),

    optparse::make_option(opt_str = c("--initial_element_mean"), action = "store", type = "double", default = 1.5,
    help = "This argument takes a floating point number (double) which defines the mean, around which the initial length/width/length measurements
    should be generated.
    By default the script sets the initial element mean to 1.5."),

    optparse::make_option(opt_str = c("--initial_element_sd"), action = "store", type = "double", default = 0.3,
    help = "This argument takes a floating point number (double) which defines the standard deviation (sd), around  the mean during initial element generation.
    By default the script sets the initial element mean to 0.3."),

    optparse::make_option(opt_str = c("--growth_variation_min"), action = "store", type = "double", default = 0.2,
    help = "This argument takes a floating point number (double) which defines the minimum growth percentage value by which each dimension measurements grow,
    from measurement to measurement. (Note: the growth value will randomly fall between the minimum and maximum values)
    By default the script sets the minimum growth variation to 0.2 (20%)."),

    optparse::make_option(opt_str = c("--growth_variation_max"), action = "store", type = "double", default = 0.5,
    help = "This argument takes a floating point number (double) which defines the maximum growth percentage value by which each dimension measurements grow,
    from measurement to measurement. (Note: the growth value will randomly fall between the minimum and maximum values)
    By default the script sets the minimum growth variation to 0.5 (50%)."),

    optparse::make_option(opt_str = c("--treatment_requested"), action = "store_true", default = FALSE, dest = "treatment",
    help = "This argument controls if the program should mimic the use of a treatment by modifying the generated measurement dimensions by a random amount,
    between a set range. \n
    By default the option is set to FALSE."),

    optparse::make_option(opt_str = c("--treatment_start"), action = "store", type = "integer", default = 3,
    help = "This argument takes an integer which defines at which measurements the treatment effect should start to be added.
    By default the script sets the number of start point to measurement 3."),

    optparse::make_option(opt_str = c("--treatment_variation_min"), action = "store", type = "double", default = -0.2,
    help = "This argument takes a floating point number (double) which defines the minimum treatment percentage value by which each dimension measurements shrinks,
    from measurement to measurement. It is important, that the treatment_variation values should be negative or close to 0 tom mimic growth arrest or shrinkage.
    (Note: the shrink value will randomly fall between the minimum and maximum values) \n
    By default the script sets the minimum treatment variation to -0.2 (-20%)."),

    optparse::make_option(opt_str = c("--treatment_variation_max"), action = "store", type = "double", default = 0.1,
    help = "This argument takes a floating point number (double) which defines the maximum treatment percentage value by which each dimension measurements shrinks,
    from measurement to measurement. It is important, that the treatment_variation values should be negative or close to 0 tom mimic growth arrest or shrinkage.
    (Note: the shrink value will randomly fall between the minimum and maximum values) \n
    By default the script sets the maximum treatment variation to 0.1 (10%)."),

    optparse::make_option(opt_str = c("--continuous_treatment"), action = "store_true", default = TRUE, dest = "continuous",
    help = "This argument controls if the growth limiting treatment effects should be applied continuously during random measurement generation, or only once,
    followed by a subtraction of treatment effect values form the already calculated growth values. The latter will mimic the effects of a one time treatment, 
    followed by a normal growth recovery.
    By default the option is set to TRUE, to mimic a continuous treatment effect."),

    optparse::make_option(opt_str = c("--format_data"), action = "store_true", default = TRUE, dest = "format",
    help = "This argument controls if the output data should be formatted to be compatible for an f-constant based tumor volume calculation CLI tool.
    (Note: the reason for this CLI tool is to generate synthetic test data for the aforementioned program, so this option will be set to TRUE by default)
    By default the option is set to TRUE."),

    optparse::make_option(opt_str = c("--treatment_group_id"), action = "store", type = "character", default = "Ctrl",
    help = "This argument controls the treatment group ID noted into the Treatment_group_id column if the format_data option was requested.
    By default the starting date is 'Ctrl'."),

    optparse::make_option(opt_str = c("--treatment_type"), action = "store", type = "character", default = "Control",
    help = "This argument controls the treatment type noted into the Treatment column, if the format_data option was requested.
    By default the starting date is 'Control'."),

    optparse::make_option(opt_str = c("--plot_data"), action = "store_true", default = TRUE, dest = "plot",
    help = "This argument controls if the output data should be plotted if the format_data option was requested.
    By default the option is set to TRUE."),

    optparse::make_option(opt_str = c("--save_plots"), action = "store_true", default = TRUE, dest = "save_plots",
    help = "This argument controls if the data plots should be saved as .png files if the format_data option and plot_data option were requested.
    By default the option is set to TRUE."),

    optparse::make_option(opt_str = c("--output_path"), action = "store", type = "character", default = def_output_path,
    help = "This argument takes a character string which defines the output path to the output .csv and .png files after generation.
    By default the path is set to be the local directory, where the .R script is located.")

)

# Create a program description
prog_descr <- c(paste0("This CLI (Command Line Interface) tool was originally created to generate synthetic sample data for the Tumor_volume_estimator script,",
" based on the f-constant formula described by Sápi et al. (2015). \nThis tool simulates the growth of subcutaneously injected tumor cells, \n",
"forming tumors on the flanks of the injected animals. The model mimics either exponential or logistical tumor growth and a treatment response, \n",
"which may vary depending on whether continuous or single-dose treatments are applied. \n",
"For standalone execution, use the syntax 'Rscript Synthetic_data_generator_V[option] [optional args]'. \n",
"Author: Gábor Bakos (UBE2C @ GitHub)"))

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
create_dates = function(no_of_days, start_date, by_days) {
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
tumor_data_generator = function(number_of_measurements, number_of_samples, growth_model, carrying_capacity, base_intrinsic_growth_rate,
    request_dates, dates, seed, initial_element_mean, initial_element_sd, growth_variation_min, growth_variation_max, treatment_administered, treatment_effect_start, 
    treatment_variation_min, treatment_variation_max, continuous_treatment) {
    
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

    #Initialize the carrying capacity and intrinsic growth rate if the logistic model was requested
    if (growth_model == "logistic") {
        intr_growth <- base_intrinsic_growth_rate

    }
    
    
    
    ##Set random number generation seed
    set.seed(seed = seed)

    
    ##Generate sample data based on exponential growth
    
    #This if statement controls which growth model should be used during data generation
    if (growth_model == "exponential") {
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
                            measurement[index] <- measurement[index] + measurement[index] * treatment_effect

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


    } else if (growth_model == "logistic") {
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
                    growth_rate <- intr_growth * (1 + growth_variation)
                    
                    #Create the first measurement of the dimension vector
                    measurement[1] <- initial_element

                    #Add the rest of the measurements which mimic the growth of the tumor by a random % range (the % range is given as the growth_variation min/max values)
                    measurement[index + 1] <- measurement[index] + growth_rate * measurement[index] * (1 - measurement[index] / carrying_capacity)
                
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
                            measurement[index] <- measurement[index] + growth_rate * measurement[index] * (1 - measurement[index] / carrying_capacity) * treatment_effect

                            #Create the rest of the measurements of the given dimension to mimic the treatment effect by a random % range
                            measurement[index + 1] <- measurement[index] + abs(growth_rate * measurement[index] * (1 - measurement[index] / carrying_capacity)) * treatment_effect

                        }

                    } else {
                        
                        #This is the innermost for loop responsible for modifying the previously generated measurements for a single dimension (either L, W or H)
                        for (index in seq(from = treatment_effect_iteration, to = nMeasurements, by = 1)) {
                            
                            #Declare the dynamic growth_variation variable (this needs to be reset for every measurement of a dimension, hence it is declared inside the loop)
                            treatment_effect <- runif(n = 1, min = treatment_variation_min, max = treatment_variation_max)

                            #Modify the rest of the measurements of the given dimension to mimic the treatment effect by a random % range
                            measurement[index] <- measurement[index] + growth_rate * measurement[index] * (1 - measurement[index] / carrying_capacity) * treatment_effect

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
        
    } else {
        stop("The growth model must be either 'exponential' or 'logistic'.")

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
format_tumor_data = function(data, treatment_group_id, treatment) {
    ##Declare function variables

    #Declare the output list
    output_list <- vector(mode = "list", length = 2)


    ##Start the dynamic variable declaration and transformation process
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

    #Name the output list elements
    names(output_list) <- names(data)
    

    ##Return the transformed dataframes
    return(output_list)

}


################################################# Section end #################################################






################################################# MARK: Subfunction - plot synthetic data #################################################
                                                #                                         #




##Transform the dataframes to a long format for plotting
plot_tumor_volumes = function(data, number_of_measurements, number_of_samples, dates, formatted_data = TRUE) {
    ##Declare/initialize function variables
    
    #Initialize the input dataframe (tumor volumes)
    input_dataframe <- data[["Tumor_volumes"]]
    
    #Declare the pivoted long dataframe
    if (formatted_data == TRUE) {
        long_df <- data.frame(matrix(nrow = nrow(input_dataframe) * number_of_measurements, ncol = 5))

    } else {
        long_df <- data.frame(matrix(nrow = nrow(input_dataframe) * number_of_measurements, ncol = 2))

    }
    

    #Declare the vector for the day x-axis breaks
    days <- vector(mode = "numeric", length = number_of_measurements)
    by_days <- vector(mode = "numeric", length = 1)


    ##Calculate the numbers for the x-axis breaks
    if (formatted_data == TRUE){
        by_days <- as.numeric(format(dates, format = "%d")[2]) - as.numeric(format(dates, format = "%d")[1])
        days <- seq(from = 0, length.out = number_of_measurements, by = by_days)
    } else {
        days <- as.numeric(unlist(stringr::str_extract_all(string = colnames(input_dataframe), pattern = "[0-9]$")))
    }
     
    

    ##Rename the dataframes columns for plotting (based on the convention of days post treatment)
    if (formatted_data == TRUE) {
        colnames(input_dataframe)[4:length(input_dataframe)] <- paste0("Day", days)

    } else {
        colnames(input_dataframe)[seq_along(colnames(input_dataframe))] <- paste0("Day", days)
        
        #Note: I must create a sample ID column for the plotting to work properly
        input_dataframe <- dplyr::mutate(.data = input_dataframe, Sample_ID = rownames(input_dataframe), .before = 1)

    }
    
    
    
    ##Pivot the input dataframe for visualization and adjust the days post treatment column
    
    #Pivot the dataframe
    if (formatted_data == TRUE) {
        long_df <- tidyr::pivot_longer(input_dataframe, cols = 4:ncol(input_dataframe), names_to = "Days_post_treatment", values_to = "Volumes")

    } else {
        long_df <- tidyr::pivot_longer(input_dataframe, cols = 2:ncol(input_dataframe), names_to = "Days_post_treatment", values_to = "Volumes")

    }
    
    
    #Adjust the Days_post_treatment_column
    long_df$Days_post_treatment <- as.numeric(stringr::str_remove_all(string = long_df$Days_post_treatment, pattern = "[a-zA-Z]"))
    print(long_df)

    ## Plot the resulting unified dataframes
    if (formatted_data == TRUE) {
        Volume_plot <- ggplot2::ggplot(data = long_df,
                    mapping = aes(x = Days_post_treatment, y = Volumes, fill = Mouse_ID, color = Mouse_ID)) +
                    ggplot2::geom_point(show.legend = TRUE) +
                    #ggplot2::geom_line() +
                    ggalt::geom_xspline(spline_shape = -0.4, show.legend = FALSE) +
                    ggplot2::scale_x_continuous(breaks = days, labels = as.character(days)) +
                    ggplot2::ggtitle("Projected tumor volumes") +
                    labs(x = expression("Days post treatment"), y = expression("Tumor volumes mm"^3),
                        color = expression("Mouse IDs"), fill = expression("Mouse IDs")) +
                    ggplot2::theme_classic() +
                    ggplot2::theme(plot.title = element_text(hjust = 0.5))

    } else {
        Volume_plot <- ggplot2::ggplot(data = long_df,
                    mapping = aes(x = Days_post_treatment, y = Volumes, fill = Sample_ID, color = Sample_ID)) +
                    ggplot2::geom_point(show.legend = TRUE) +
                    #ggplot2::geom_line() +
                    ggalt::geom_xspline(spline_shape = -0.4, show.legend = FALSE) +
                    ggplot2::scale_x_continuous(breaks = days, labels = as.character(days)) +
                    ggplot2::ggtitle("Projected tumor volumes") +
                    labs(x = expression("Days post treatment"), y = expression("Tumor volumes mm"^3),
                        color = expression("Sample IDs"), fill = expression("Sample IDs")) +
                    ggplot2::theme_classic() +
                    ggplot2::theme(plot.title = element_text(hjust = 0.5))
    }
    

    return(Volume_plot)
}


################################################# Section end #################################################






################################################# MARK: Main function #################################################
                                                #                     #




##The main function integrating all of the subfunctions
main = function(measurements = arguments$number_of_measurements, samples = arguments$number_of_samples, model = arguments$growth_model, carry_cap = arguments$carrying_capacity, 
    base_intr_growth_rate = arguments$base_intrinsic_growth_rate, rdates = arguments$dates, sdate = arguments$start_date, bydays = arguments$by_days, set_seed = arguments$seed,
    init_element_mean = arguments$initial_element_mean, init_element_sd = arguments$initial_element_sd, gr_variation_min = arguments$growth_variation_min, gr_variation_max = arguments$growth_variation_max,
    treatment = arguments$treatment, treatment_start = arguments$treatment_start, tr_variation_min = arguments$treatment_variation_min, tr_variation_max = arguments$treatment_variation_max,
    cont_treatment = arguments$continuous, format = arguments$format, tr_group_id = arguments$treatment_group_id, tr_type = arguments$treatment_type,
    plot = arguments$plot, save_pl = arguments$save_plots, output_path = arguments$output_path) {
        ##Call the package management subfunctions 
        package_controller()
        suppressPackageStartupMessages(package_loader())


        ##Call the create dates subfunction if requested
        if (rdates == TRUE) {
            Dates <- create_dates(no_of_days = measurements, start_date = sdate, by_days = bydays)

        }
        

        ##Call the tumor_data-generator subfunction
        message("Generating synthetic data... \n")
        
        Synthetic_data <- tumor_data_generator(number_of_measurements = measurements, number_of_samples = samples, growth_model = model, carrying_capacity = carry_cap,
            base_intrinsic_growth_rate = base_intr_growth_rate, request_dates = rdates, dates = Dates, seed = set_seed, initial_element_mean = init_element_mean,
            initial_element_sd = init_element_sd, growth_variation_min = gr_variation_min, growth_variation_max = gr_variation_max, treatment_administered = treatment,
            treatment_effect_start = treatment_start, treatment_variation_min = tr_variation_min, treatment_variation_max = tr_variation_max, continuous_treatment = cont_treatment)

        
        ##Call the format data subfunction if requested
        if (format == TRUE) {
            message("Data formatting was requested. \n")

            Formatted_synth_data <- format_tumor_data(Synthetic_data, treatment_group_id = tr_group_id, treatment = tr_type)

        }


        ##Call the plot_tumor_volumes subfunction if requested
        if (plot == TRUE && format == TRUE) {
            message("Data plotting was requested. Generating data plot... \n")

            Synth_data_plot <- plot_tumor_volumes(data = Formatted_synth_data, number_of_measurements = measurements, number_of_samples = samples, dates = Dates)
            
        } else {
            message("Data plotting was requested. Generating data plot... \n")

            Synth_data_plot <- plot_tumor_volumes(data = Synthetic_data, number_of_measurements = measurements, number_of_samples = samples, dates = Dates)

        }


        ##Save the plots if requested
            message("Saving the data plot... \n")
                
            if (save_pl == TRUE) {
                #NOTE: I'm using an older save syntax and I'm sending the file directly to the device as the ggsave does not seem to work!
                png(filename = paste0("Synthetic_tumor_volumes_plot.png"), width = 1500, height = 750, units = "px")
                print(Synth_data_plot)
                dev.off()
            }

        
        ##Save the generated synthetic data or formatted data if requested
        message("Saving the generated synthetic data... \n")

        

        if (format == TRUE) {
            for (element in seq_along(Formatted_synth_data)) {
                readr::write_csv(x = Formatted_synth_data[[element]], file = paste0(output_path, "/", "Form_", names(Formatted_synth_data)[[element]], ".csv"))
            }
            
        } else {
            for (element in seq_along(Synthetic_data)) {
               readr::write_csv(x = Synthetic_data[[element]], file = paste0(output_path, "/", names(Synthetic_data)[[element]], ".csv"))
            }

        }

        
        

        
        ##Success message
        message("The data generation was successful. \n",
        "The generated .csv and .png files were saved to the ", arguments$output_path, " directory. \n")


        ##Print session inf
        print(sessionInfo())
}


################################################# Section end #################################################






################################################# MARK: Main function call #################################################
                                                #                          #




## Call the main function as an execution point for the script
main()


################################################# Section end #################################################