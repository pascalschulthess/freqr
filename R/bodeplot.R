#' Bode plot
#'
#' @param fr Frequency range in powers of 10
#' @param n Numerator coefficients of a transfer function
#' @param d Denominator coefficients of a transfer function
#' @param A Matrix A of a state-space representation
#' @param B Matrix B of a state-space representation
#' @param C Matrix C of a state-space representation
#' @param D Matrix D of a state-space representation
#' @param highlight Optional highlights in the plot
#' @param label Optional labels on the highlights
#' @return Bode magnitude and phase plot
#' @examples
#' bodeplot(c(0, 1), n = c(1, 0.1, 7.5), d = c(1, 0.12, 9, 0, 0))
#' bodeplot(c(0, 1), A = matrix(c(-1,1,0,-1), nrow = 2),
#'                   B = matrix(c(-1,1), nrow = 2),
#'                   C = matrix(c(0,1), ncol=2),
#'                   D = 0, highlight = 2, label = "Highlight")
#'
#' @export

bodeplot <- function(fr, n = NULL, d = NULL, A = NULL, B = NULL, C = NULL, D = NULL, highlight = NULL, label = NULL) {

  # Check inputs ------------------------------------------------------------
  # Check method
  if (!is.null(n) && !is.null(d)) {
    method = "tf"
  } else if (!is.null(A) && !is.null(B) && !is.null(C) && !is.null(D)) {
    method = "state-space"
  }

  # Check numerator and denominator
  if (method == "tf") {
    if (length(n) == 0) {
      stop("Numerator is empty")
    } else if (length(d) == 0) {
      stop("Denominator is empty")
    }
  }

  # Check matrices A, B, C, and D
  if (method == "state-space") {
    if (is.numeric(A) == FALSE | is.numeric(A) == FALSE | is.numeric(A) == FALSE | is.numeric(A) == FALSE) {
      stop("Matrices A, B, C, and D must be numeric")
    }
    if (dim(A)[1] != dim(A)[2]) {
      stop("Matrix A must be square")
    } else if (dim(B)[2] != 1) {
      stop("Matrix B must have dimension: n x 1. For MIMO model use the sigmaplot function")
    } else if (dim(C)[1] != 1) {
      stop("Matrix C must have dimension: 1 x n. For MIMO model use the sigmaplot function")
    } else if (!is.null(dim(D))) {
      stop("Matrix D must have dimension: 1 x 1. For MIMO model use the sigmaplot function")
    }
  }

  # Check labels
  if (!is.null(highlight) && !is.null(label) && label[1] != "numeric" && length(highlight) != length(label)) {
    stop("Custom labels must have the same length as the highlights.")
  }


  # Frequency ---------------------------------------------------------------
  f <- 10 ^ seq(from = fr[1], to = fr[2], length = 1000)


  # Calculate transfer function ---------------------------------------------
  # Coefficients of numerator and denominator
  if (method == "tf") {

    # Transfer function
    G <- pracma::polyval(n, complex(imaginary = f))/pracma::polyval(d, complex(imaginary = f))
  }

  # A, B, C, D matrices of the linear state space
  if (method == "state-space") {

    # Initialise G
    G <- rep(NA, length(f))

    # Calculate G from state space matrices
    for (i in 1:length(f)) {
      G[i] <- C %*% solve(complex(imaginary = f[i])*diag(dim(A)[1]) - A) %*% B + D
    }
  }


  # Calculate magnitude and phase -------------------------------------------
  # Magnitude
  mag <-  20*log10(Mod(G))

  # Phase
  phase <- 180/pi*Arg(G)

  # If there's a sign switch in the phase
  signswitch <- which(diff(sign(phase)) != 0)
  if (length(signswitch) != 0) {
    if (sign(phase[signswitch]) == 1) {
      phase[1:signswitch] <- phase[1:signswitch] - 360
    } else if (sign(phase[signswitch]) == -1) {
      phase[1:signswitch] <- phase[1:signswitch] + 360
    }
  }

  # Put everything in a dataframe
  res <- data.frame(Frequency = f,
                    Magnitude = mag,
                    Phase = phase)

  # Add labels if desired ---------------------------------------------------
  if (!is.null(highlight)) {
    if (length(label) == 1 && label[1] == "numeric") {
      label.mag <- paste(round(log10(f), 2), round(mag, 2), sep = ", ")
      label.mag <- paste0("10^", label.mag)
      label.phase <- paste(round(log10(f), 2), round(phase, 2), sep = ", ")
      label.phase <- paste0("10^", label.phase)

      res$LabelMag = as.character(label.mag)
      res$LabelPhase = as.character(label.phase)
    } else if (!is.null(label) && label != "numeric") {
      lab <- rep(NA, length(f))
      for (i in 1:length(highlight)) {
        h <- which.min(abs(res$Frequency - highlight[i]))
        lab[h] <- label[i]
      }
      res$LabelMag = as.character(lab)
      res$LabelPhase = as.character(lab)
    }
  }

  # Plot --------------------------------------------------------------------
  # Magnitude
  p1 <- ggplot2::ggplot(data = res, ggplot2::aes_string(x = "Frequency", y = "Magnitude")) +
    ggplot2::geom_line(ggplot2::aes_string(x = "Frequency", y = "Magnitude")) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed") +
    ggplot2::scale_x_log10(breaks = 10 ^ seq(fr[1], fr[2], 1)) +
    ggplot2::ggtitle("Bode diagram") +
    ggplot2::ylab("Magnitude [dB]")

  # Phase
  p2 <- ggplot2::ggplot(data = res, ggplot2::aes_string(x = "Frequency", y = "Phase")) +
    ggplot2::geom_line(ggplot2::aes_string(x = "Frequency", y = "Phase")) +
    ggplot2::scale_y_continuous(breaks = seq(from = -360, to = 360, by = 45)) +
    ggplot2::scale_x_log10(breaks = 10 ^ seq(fr[1], fr[2], 1), labels = scales::comma) +
    ggplot2::xlab('Frequency [rad/s]') +
    ggplot2::ylab("Phase [deg]")

  # Add highlights and labels if desired ------------------------------------
  if (!is.null(highlight)) {

    # Highlights follow the ggplot2 color scheme
    hues = seq(15, 375, length = length(highlight) + 1)
    cols <- grDevices::hcl(h = hues, l = 65, c = 100)[1:length(highlight)]

    # Add highlights to plot
    for (i in 1:length(highlight)) {
      h <- which.min(abs(res$Frequency - highlight[i]))
      p1 <- p1 + ggplot2::geom_point(data = res[h,, drop = FALSE], mapping = ggplot2::aes_string(x = "Frequency", y = "Magnitude"), color = cols[i], size = 3)
      p2 <- p2 + ggplot2::geom_point(data = res[h,, drop = FALSE], mapping = ggplot2::aes_string(x = "Frequency", y = "Phase"), color = cols[i], size = 3)

      if (!is.null(label)) {
        p1 <- p1 + ggrepel::geom_label_repel(data = res[h,, drop = FALSE], mapping = ggplot2::aes_string(x = "Frequency", y = "Magnitude", label = "LabelMag"), fill = cols[i], color = "white", box.padding = grid::unit(0.25, "lines"),
                                             point.padding = grid::unit(0.5, "lines"))
        p2 <- p2 + ggrepel::geom_label_repel(data = res[h,, drop = FALSE], mapping = ggplot2::aes_string(x = "Frequency", y = "Phase", label = "LabelPhase"), fill = cols[i], color = "white", box.padding = grid::unit(0.25, "lines"),
                                             point.padding = grid::unit(0.5, "lines"))
      }
    }
  }

  # Combine everything
  grid::grid.newpage()
  grid::grid.draw(rbind(ggplot2::ggplotGrob(p1), ggplot2::ggplotGrob(p2), size = 'last'))
}
