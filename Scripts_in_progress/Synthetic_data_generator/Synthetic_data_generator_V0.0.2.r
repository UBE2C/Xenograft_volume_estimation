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

measurement_generator = function(dimension_length, growth_variation_min = 0.2, growth_variation_max = 0.5, initial_element_mean = 1, initial_element_sd = 1) {
    initial_element <- vector(mode = "numeric", length = 1)
    growth_variation <- vector(mode = "numeric", length = 1)
    vector_length <- dimension_length
    dimension <- vector(mode = "numeric", length = dimension_length)

    while (initial_element <= 0.5) {
        initial_element <- rnorm(n = 1, mean = initial_element_mean, sd = initial_element_sd)
    }


    for (index in seq_len(vector_length - 1)){
        growth_variation <- runif(n = 1, min = growth_variation_min, max = growth_variation_max)
        
        dimension[1] <- initial_element
        dimension[index + 1] <- dimension[index] + dimension[index] * growth_variation
        
    }

    return(dimension)
}


sample_generator = function(number_of_measurements) {
    dimensions <- c("L", "W", "H")
    sample <- list()

    for (dim in dimensions) {
        sample[[dim]] <- measurement_generator(number_of_measurements)
    }

    return(sample)
}





sample_generator = function(number_of_measurements, growth_variation_min = 0.2, growth_variation_max = 0.5, initial_element_mean = 1, initial_element_sd = 1) {
    
    growth_variation <- vector(mode = "numeric", length = 1)
    vector_length <- number_of_measurements
    measurement <- vector(mode = "numeric", length = number_of_measurements)
    dimensions <- c("L", "W", "H")
    sample <- list()

    


    

    for (dim in dimensions) {
        initial_element <- vector(mode = "numeric", length = 1)
        
        while (initial_element <= 0.5) {
        initial_element <- rnorm(n = 1, mean = initial_element_mean, sd = initial_element_sd)
        }
        
        for (index in seq_len(vector_length - 1)){
        growth_variation <- runif(n = 1, min = growth_variation_min, max = growth_variation_max)
        
        measurement[1] <- initial_element
        measurement[index + 1] <- measurement[index] + measurement[index] * growth_variation
        
        }

        sample[[dim]] <- measurement
    }

    return(sample)
}



## Calculate tumor volumes for a sample
volume_data_generator = function(sample) {
    #sample_volumes <- list()
    number_of_measurements <- length(sample[[1]])
    volume <- vector(mode = "numeric", length = number_of_measurements)

    for (index in seq_len(number_of_measurements)) {
        volume[index] <- (4/3) * pi * (sample$L[index]/2) * (sample$W[index]/2) * (sample$H[index]/2)

    }

    return(volume)    

}



tumor_data_generator = function(number_of_measurements, number_of_samples, dates, seed = 1234,
    growth_variation_min = 0.2, growth_variation_max = 0.5, initial_element_mean = 1, initial_element_sd = 1) {
    
    ##Declare function variables

    #A vector containing the factor for tumor growth variation
    growth_variation <- vector(mode = "numeric", length = 1)

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








