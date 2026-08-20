#' Get first and total Sobol indices by Carlos Camino
#'
#' This function calculates the first-order and total Sobol indices for sensitivity analysis
#' based on the input data and predictions.
#'
#' @param data A data frame containing the input variables and the predicted output.
#' @param output A character string specifying the name of the column with the predictions.
#' @param N An integer specifying the number of samples to be used for Sobol sensitivity analysis.
#' @param normalize A boolean value; if TRUE, the Sobol indices will be normalized.
#'
#' @return A data frame containing the first-order and total Sobol indices for each input variable.
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage
#' data <- data.frame(input1 = runif(100), input2 = runif(100), output = rnorm(100))
#' sobol_indices <- get.sobol.indices(data, output = "output", N = 1000, normalize = TRUE)
#' print(sobol_indices)
#' }
get.sobol.indices <- function(data, output, N,normalize = FALSE) {
  # Extract the input variables
  inputs_ <- setdiff(names(data), output)

  if (N > nrow(data) / 2) {
    stop("N cannot be greater than or equal to half the number of rows in the dataset.")
  }

  # Initialize a matrix to store the Sobol indices
  sobol.1 <- matrix(0, nrow = 2, ncol = length(inputs_),
                         dimnames = list(c("Si", "STi"), inputs_))
  sobol.2 <-sobol.1

  # Normalize the input variables
  normalized_data <- scale(data[, inputs_])
  # Combine the normalized input variables with the output variable
  normalized_data <- cbind(normalized_data, data[, output])
  colnames(normalized_data)<- c(inputs_,output)
  rownames(normalized_data) <- rownames(data)
  # Randomly sample N rows for the first dataset
  m.1 <- normalized_data[sample(nrow(normalized_data), N, replace = FALSE), ]

  # Select the remaining rows for the second dataset
  m.2 <- normalized_data[setdiff(1:nrow(normalized_data), rownames(m.1)), ]
  m.2<- m.2[1:N,]
  # Computation of Johnson index
  # ----------------------------------------------------------------
  ## create data.frame for Johson index
  df.johnson<- data.frame(Band = character(), Parameter=character(), Index = numeric())

  ind.johnson <- sensitivity:: johnson(normalized_data[,inputs_], y=normalized_data[,output], logistic = F) #, nboot = 100, conf=0.95)
  df.johnson <- rbind(df.johnson, data.frame(Band = output, Parameter =rownames(ind.johnson$johnson), I.Johnson = ind.johnson$johnson$original ))


  # Computation of Sobo for two dataset and merge by mean
  # ----------------------------------------------------------------

  for (i in 1:length(inputs_)) {
    # Define the model inputs (excluding the current input)
    inputs.i <- inputs_[i]


    # Computation of E(Y), V(Y), Si and Ti (dataset1)
    # ----------------------------------------------------------------

    f0 <- (1 / N) * sum(m.1[, inputs.i] * m.1[, output])
    VY <- 1 / (2 * N - 1) * sum(m.1[, inputs.i]^2 + m.1[, output]^2) - f0
    #first-order sensitivity index (Si)
    Si <- (1 / (N - 1) * sum(m.1[, inputs.i] * m.1[, output]) - f0) / VY
    #Total-order sensitivity index (STi)
    STi <- 1 - (1 / (N - 1) * sum(m.1[, output] * m.1[, output]) - f0) / VY
    sobol.1[, inputs.i] <- c(Si, STi)

    # Computation of E(Y), V(Y), Si and Ti  (dataset2)
    # ----------------------------------------------------------------

    f0 <- (1 / N) * sum(m.2[, inputs.i] * m.2[, output])
    VY <- 1 / (2 * N - 1) * sum(m.2[, inputs.i]^2 + m.2[, output]^2) - f0
    #first-order sensitivity index (Si)
    Si <- (1 / (N - 1) * sum(m.2[, inputs.i] * m.2[, output]) - f0) / VY
    #Total-order sensitivity index (STi)
    STi <- 1 - (1 / (N - 1) * sum(m.2[, output] * m.2[, output]) - f0) / VY
    sobol.2[, inputs.i] <- c(Si, STi)
  } # end for

  # Compute the average Sobol indices across both datasets
  average_sobol_indices <- (sobol.1 + sobol.2) / 2

  # Create a data frame for average Sobol indices
  average_sobol_df <- data.frame(Band = rep(output, length(inputs_)),
                                 Parameter = inputs_,
                                 Si = (average_sobol_indices[1, ]),
                                 STi = average_sobol_indices[2, ],
                                 row.names = NULL)

  # Merge the data frames based on 'Band' and 'Parameter'
  merged_df <- dplyr::left_join(average_sobol_df, df.johnson, by = c("Band", "Parameter")) |>
    dplyr::mutate(Difference = Si - I.Johnson)

  # Calculate the threshold based on the 10th percentile of the absolute values of the Johnson indices
  threshold <- quantile(abs(merged_df$I.Johnson), probs = 0.25)
  # Apply condition to set Sobol indices to zero where Johnson index is very low
  merged_df <- merged_df |>
    dplyr::mutate(Si = ifelse(abs(I.Johnson) < threshold, 0, Si),
           STi = ifelse(abs(I.Johnson) < threshold, 0, STi))

  if (normalize == TRUE) {
    merged_df <- merged_df |>
      dplyr::mutate(Si_norm = abs(Si) / sum(abs(Si)) * 100)  |>
      dplyr::mutate(I.Johnson_norm = abs(I.Johnson) / sum(abs(I.Johnson)) * 100)
      #dplyr::mutate(Si_norm = ifelse(Si != 0, abs(Si) / sum(abs(Si)) * 100, 0),
       #      I.Johnson_norm = abs(I.Johnson) / sum(abs(I.Johnson)) * 100)

    return(merged_df)
  } else {
    return(merged_df)
  }


}

