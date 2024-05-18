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



volume_data_generator = function(number_of_measurements, number_of_samples, 
    growth_variation_min = 0.2, growth_variation_max = 0.5, initial_element_mean = 1, initial_element_sd = 1) {
    
    ##Declare function variables
    growth_variation <- vector(mode = "numeric", length = 1)
    nMeasurements <- vector(mode = "numeric", length = 1)
    dimensions <- vector(mode = "character", length = 3)
    measurement <- vector(mode = "numeric", length = number_of_measurements)
    nSamples <- vector(mode = "numeric", length = 1)
    sample <- list()
    sample_volumes <- list()
    
    
    #Initialize function variables
    nMeasurements <- number_of_measurements
    nSamples <- number_of_samples
    dimensions <- c("L", "W", "H")
    

    


    for (element in seq_len(nSamples)) {
        volume <- vector(mode = "numeric", length = nMeasurements)

        for (dim in dimensions) {
            initial_element <- vector(mode = "numeric", length = 1)
        
            while (initial_element <= 0.5) {
            initial_element <- rnorm(n = 1, mean = initial_element_mean, sd = initial_element_sd)
            }
        
            for (index in seq_len(nMeasurements - 1)){
                growth_variation <- runif(n = 1, min = growth_variation_min, max = growth_variation_max)
                
                measurement[1] <- initial_element
                measurement[index + 1] <- measurement[index] + measurement[index] * growth_variation
            
            }

            sample[[dim]] <- measurement
        }

        for (ind in seq_len(nMeasurements)) {
            volume[ind] <- (4/3) * pi * (sample$L[ind]/2) * (sample$W[ind]/2) * (sample$H[ind]/2)

        }

        sample_volumes[[element]] <- volume

    }    

    return(sample_volumes)
}








