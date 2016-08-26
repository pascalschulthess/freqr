#' Sigma plot
#'
#' @param fr Frequency range in powers of 10
#' @param fr.unit Unit of the frequency
#' @param A Matrix A of a state-space representation
#' @param B Matrix B of a state-space representation
#' @param C Matrix C of a state-space representation
#' @param D Matrix D of a state-space representation
#' @param as.period Flag to plot as period
#' @return Singular value plot
#' @examples
#' sigmaplot(fr = c(-2, 2),
#'            A = matrix(c(0, 1, -1, -1), nrow = 2),
#'            B = diag(2),
#'            C = matrix(c(10, 0, 0, 1), nrow = 2),
#'            D = 0.1*diag(2))
#'
#' @export

sigmaplot <- function(fr, fr.unit = NULL,
                      A = NULL, B = NULL, C = NULL, D = NULL,
                      as.period = F) {

  # Check inputs ------------------------------------------------------------
  # Check frequency unit
  if (!is.null(fr.unit) && !is.character(fr.unit)) {
    stop("fr.unit must be a string")
  }

  # Check matrices A, B, C, and D
  if (is.numeric(A) == FALSE | is.numeric(A) == FALSE | is.numeric(A) == FALSE | is.numeric(A) == FALSE) {
    stop("Matrices A, B, C, and D must be numeric")
  }
  if (dim(A)[1] != dim(A)[2]) {
    stop("A must be square")
  } else if (dim(B)[1] != dim(A)[2]) {
    stop("When dim(A) = n x n, then dim(B) = n x p")
  } else if (dim(C)[2] != dim(A)[1]) {
    stop("When dim(A) = n x n, then dim(C) = q x n")
  } else if (!is.null(dim(D))) {
    if (all(dim(D) != c(dim(C)[1], dim(B)[2]))) {
      stop("When dim(B) = n x p and dim(C) = q x n, then dim(D) = q x p")
    }
  }

  # Check as.period
  if (!is.logical(as.period)) {
    stop("as.period must be logical")
  }

  # Frequency ---------------------------------------------------------------
  f <- 10 ^ seq(from = fr[1], to = fr[2], length = 1000)

  # Calculate singular values -----------------------------------------------
  # Initialise singular values
  sigma.max <- rep(NA, length(f))
  sigma.min <- rep(NA, length(f))

  # Calculate G and sigmas from state space matrices
  for (i in 1:length(f)) {
    G <- C %*% solve(complex(imaginary = f[i])*diag(dim(A)[1]) - A) %*% B + D
    sigmas <- svd(G)$d
    sigma.max[i] <- 20*log10(Mod(utils::head(sigmas, 1)))
    sigma.min[i] <- 20*log10(Mod(utils::tail(sigmas, 1)))
  }

  # Put everything in a dataframe
  if (!isTRUE(as.period)) {
    res <- data.frame(Frequency = f,
                      Sigma.max = sigma.max,
                      sigma.min = sigma.min)
  } else {
    res <- data.frame(Frequency = 1/f,
                      Sigma.max = sigma.max,
                      sigma.min = sigma.min)
  }

  # Plot --------------------------------------------------------------------
  # Frequency units
  if (is.null(fr.unit) && !isTRUE(as.period)) {
    fr.unit <- "rad/s"
  } else if (is.null(fr.unit) && isTRUE(as.period)) {
    fr.unit <- "s"
  }

  # Reshape data for ggplot2
  res.plt <- reshape2::melt(res,  id.vars = "Frequency", variable.name = "Sigma")

  # Magnitude plot
  p <- ggplot2::ggplot(data = res.plt, ggplot2::aes_string(x = "Frequency", y = "value")) +
    ggplot2::geom_line(ggplot2::aes_string(colour = "Sigma")) +
    ggplot2::scale_color_discrete(labels = c("max(sigma)", "min(sigma)")) +
    ggplot2::scale_x_log10(breaks = 10 ^ seq(fr[1], fr[2], 1)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed") +
    ggplot2::ggtitle("Sigma plot") +
    ggplot2::xlab("Frequency [rad/s]") +
    ggplot2::ylab("Magnitude [dB]")

  # x-axis scaling and label
  if (!isTRUE(as.period)) {
    p <- p + ggplot2::scale_x_log10(breaks = 10 ^ seq(fr[1], fr[2], 1)) +
      ggplot2::xlab(paste("Frequency", paste0("[", fr.unit, "]")))
  } else {
    p <- p + ggplot2::scale_x_log10() +
      ggplot2::xlab(paste("Period", paste0("[", fr.unit, "]")))
  }

  # Print
  print(p)
}
