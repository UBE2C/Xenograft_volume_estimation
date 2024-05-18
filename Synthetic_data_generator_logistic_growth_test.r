#!usr/bin/env -S Rscript --vanilla



#########################################
#                                       #
#   Tumor measurement data processing   #
#       Synthetic data generation       #
#                                       #
#########################################






################################################# MARK: Build partial functions #################################################
                                                #                               #




## Check for the required packages and install them if missing

dates <- c("01.01.2024", "05.01.2024", "09.01.2024", "13.01.2024", "17.01.2024", "21.01.2024")

tumor_data_generator = function(number_of_measurements, number_of_samples, dates, seed = 1234, initial_element_sd = 0.3,
    growth_variation_min = 0.2, growth_variation_max = 0.5, initial_element_mean = 1.5, treatment_variation_min = -0.1, treatment_variation_max = 0.1,
    treatment_effect_start = 3, treatment_administered = FALSE, continuous_treatment = TRUE) {
    
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
    measurement_dates <- vector(mode = "character", length = length(dates))

    #The output dataframes containing the generated synthetic tumor volume or LxW measurements
    output_volume_dataframe <- data.frame(matrix(nrow = number_of_samples, ncol = length(dates)))
    output_LxW_dataframe <- data.frame(matrix(nrow = number_of_samples, ncol = length(dates)))

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

    #Initialize the measurement dates variable
    measurement_dates <- dates

    carrying_capacity <- 15
    
    
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
                growth_rate <- 0.1 * (1 + growth_variation)
                
                #Create the first measurement of the dimension vector
                measurement[1] <- initial_element

                #Add the rest of the measurements which mimic the growth of the tumor by a random % range (the % range is given as the growth_variation min/max values)
                measurement[index + 1] <- measurement[index] + growth_rate * measurement[index] * (1 - measurement[index]/carrying_capacity)
            
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
                        measurement[index] <- measurement[index] + growth_rate * measurement[index] * (1 - measurement[index]/carrying_capacity) * treatment_effect

                        #Create the rest of the measurements of the given dimension to mimic the treatment effect by a random % range
                        measurement[index + 1] <- measurement[index] + abs(growth_rate * measurement[index] * (1 - measurement[index]/carrying_capacity)) * treatment_effect

                    }

                } else {
                    
                    #This is the innermost for loop responsible for modifying the previously generated measurements for a single dimension (either L, W or H)
                    for (index in seq(from = treatment_effect_iteration, to = nMeasurements, by = 1)) {
                        
                        #Declare the dynamic growth_variation variable (this needs to be reset for every measurement of a dimension, hence it is declared inside the loop)
                        treatment_effect <- runif(n = 1, min = treatment_variation_min, max = treatment_variation_max)

                        #Modify the rest of the measurements of the given dimension to mimic the treatment effect by a random % range
                        measurement[index] <- measurement[index] + growth_rate * measurement[index] * (1 - measurement[index]/carrying_capacity) * treatment_effect

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

    #Name the rows and columns of the output dataframes
    colnames(output_volume_dataframe) <- measurement_dates
    rownames(output_volume_dataframe) <- sample_names
    colnames(output_LxW_dataframe) <- measurement_dates
    rownames(output_LxW_dataframe) <- sample_names

    #Add the output dataframe to the final output list
    final_output_list[[1]] <- output_volume_dataframe
    final_output_list[[2]] <- output_LxW_dataframe

    #Name the final output list elements
    names(final_output_list) <- c("Tumor_volumes", "Tumor_LxW_values")
    

    ##Return the output dataframe
    return(final_output_list)

}








