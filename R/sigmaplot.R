#' Sigma plot
#'
#' @param fr Frequency range in powers of 10
#' @param A Matrix A of a state-space representation
#' @param B Matrix B of a state-space representation
#' @param C Matrix C of a state-space representation
#' @param D Matrix D of a state-space representation
#' @return Singular value plot
#' @examples
#' sigmaplot(c(-2, 2), A <- matrix(c(0, 1, -1, -1), nrow = 2)
#'                     B <- diag(2)
#'                     C <- matrix(c(10, 0, 0, 1), nrow = 2)
#'                     D <- 0.1*diag(2))

sigmaplot <- function(fr,
                      A = NULL, B = NULL, C = NULL, D = NULL) {

  # Frequency
  f <- 10 ^ seq(from = fr[1], to = fr[2], length = 1000)

  # Initialise singular values
  sigma.max <- rep(NA, length(f))
  sigma.min <- rep(NA, length(f))

  # Calculate G and sigmas from state space matrices
  for (i in 1:length(f)) {
    G <- C %*% solve(complex(imaginary = f[i])*diag(dim(A)[1]) - A) %*% B + D
    sigmas <- svd(G)$d
    sigma.max[i] <- head(sigmas, 1)
    sigma.min[i] <- tail(sigmas, 1)
  }

  # Plotting preparations (sigmas in dB)
  df <- data.frame(frequency = f,
                   sigma.max = 20*log10(Mod(sigma.max)),
                   sigma.min = 20*log10(Mod(sigma.min)))

  res <- reshape2::melt(df,  id.vars = "frequency", variable.name = "sigma")

  # Magnitude plot
  ggplot2::ggplot(res, ggplot2::aes(frequency, value)) +
    ggplot2::geom_line(ggplot2::aes(colour = sigma)) +
    ggplot2::scale_color_discrete(labels = c("max(sigma)", "min(sigma)")) +
    ggplot2::scale_x_log10(breaks = 10 ^ seq(fr[1], fr[2], 1)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed") +
    ggplot2::ggtitle("Sigma plot") +
    ggplot2::xlab("Frequency [rad/s]") +
    ggplot2::ylab("Magnitude [dB]") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
          legend.position = "bottom")
}