update_scHOT_output <- function(scHOT, new_output, output_name) {

  scHOT_output <- scHOT@scHOT_output


  # if previous is NULL

  if (nrow(scHOT_output) == 0 | ncol(scHOT_output) == 0) {

    scHOT_output <- DataFrame(new_output,
                              row.names = paste(scHOT@testingScaffold[, 1],
                                               scHOT@testingScaffold[, 2], sep = "_"))
    colnames(scHOT_output) <- output_name

  } else {

    current_testing <-  paste(scHOT@testingScaffold[, 1], scHOT@testingScaffold[, 2], sep = "_")

    # if (class(new_output) == "list") {
    #   scHOT_output[current_testing, ][[output_name]] <- new_output
    # } else {
    #   scHOT_output[current_testing, output_name] <- new_output
    # }
    scHOT_output[current_testing, output_name] <- new_output

  }


  scHOT@scHOT_output <- scHOT_output

  return(scHOT)
}



update_stratifiedSample <- function(scHOT, numberPermutations, numberScaffold) {

  current_testing <-  paste(scHOT@testingScaffold[, 1], scHOT@testingScaffold[, 2], sep = "_")

  scHOT@scHOT_output[current_testing, ]$numberPermutations[stratifiedSample(
    scHOT@scHOT_output[current_testing, ]$globalHigherOrderFunction,
    length = numberScaffold)] <- numberPermutations

  return(scHOT)

}
