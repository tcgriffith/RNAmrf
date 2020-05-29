#' Title gg_mat
#'
#' Use ggplot to visualize matrix
#'
#' @param mat a matrix
#'
#' @return ggplot2 object
#' @export
#' @import ggplot2
#'
gg_mat = function(mat) {
  df = as.data.frame(which(!is.na(mat), arr.ind = TRUE))
  df$val = mat[which(!is.na(mat))]
  ggplot(df) + geom_tile(aes(x=col, y=row, fill = val)) + coord_fixed()
}
