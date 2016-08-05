#' Nyquist plot
#'
#' @param fr Frequency range in powers of 10
#' @param n Numerator coefficients of a transfer function
#' @param d Denominator coefficients of a transfer function
#' @return Nyquist plot
#' @examples
#' nyquistplot(c(-2, 3), n = c(2, 5, 1), d = c(1, 2, 3))

nyquistplot <- function(fr, n = NULL, d = NULL) {

  # Frequency
  f <- 10 ^ seq(from = fr[1], to = fr[2], length = 1000)

  # Transfer function
  tf <- pracma::polyval(n, complex(imaginary = f))/pracma::polyval(d, complex(imaginary = f))

  # Unit circle
  unit <- complex(real = cos(seq(0, 2*pi, length = 1000)), imaginary = sin(seq(0, 2*pi, length = 1000)))

  # Make dataframe
  res <- data.frame(TransferFunction = tf,
                   Unit = unit)

  # Plot
  ggplot2::ggplot(res) +
    ggplot2::geom_path(ggplot2::aes(Re(Unit), Im(Unit)), color = "grey") +
    ggplot2::geom_point(ggplot2::aes(-1,0), color = "red") +
    ggplot2::geom_path(ggplot2::aes(Re(TransferFunction), Im(TransferFunction)), arrow = grid::arrow(angle = 15, ends = "last", type = "closed")) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::ggtitle("Nyquist plot") +
    ggplot2::xlab("Real Axis") +
    ggplot2::ylab("Imaginary Axis") +
    ggplot2::theme_bw()
}
