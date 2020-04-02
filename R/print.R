#' @export
print.mediator <- function(x, digits = 5, ...){
  clean_x <- dplyr::as_tibble(x) %>%
    dplyr::mutate_if(is.numeric, ~round(., digits))

  print(clean_x)
}
