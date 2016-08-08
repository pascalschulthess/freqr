# freqr: Frequency-domain response analysis in R

## Summary

A collection of `R` functions that enables the frequency-domain analysis of dynamic systems.

## Installation

```
devtools::install_github("pascalschulthess/freqr")
```

## Functions and their usage

### Bode plot

Takes the coefficients of a transfer function or the matrices A, B, C and D of a linear SISO state space representation and generates a Bode plot.

#### Usage

```
bodeplot <- function(fr, n = NULL, d = NULL, 
                 A = NULL, B = NULL, C = NULL, D = NULL,
                 highlight = NULL, label = NULL)
```

#### Example

In order to plot the transfer function

```
		 s^2 + 0.1s + 7.5
G(s) = ----------------------
		s^4 + 0.12s^3 + 9s^2
```

from 10^0 to 10^1 and add a highlight at 2 [rad/s], the function inputs take the form

```
bodeplot(fr = c(0, 1),
	      n = c(1, 0.1, 7.5),
	      d = c(1, 0.12, 9, 0, 0), 
	      highlight = 2,
	      label = "Highlight")
```

and result in

![](https://github.com/pascalschulthess/freqr/blob/master/readme-examples/bodeplot-example.png)

#### Dependencies

```
library(pracma)
library(ggplot2)
library(dplyr)
library(grid)
```

### Nyquist plot

Takes the coefficients of a transfer function and generates a Nyquist plot.

#### Usage

```
nyquistplot(fr, n = NULL, d = NULL)
```

#### Example

In order to plot the transfer function

```
		     100
G(s) = ----------------
		s^2 + 6s + 100
```

from 10^-2 to 10^3 the function inputs take the form

```
nyquistplot(fr = c(-2, 3), n = 100, d = c(1, 6, 100))
```

and result in

![](https://github.com/pascalschulthess/freqr/blob/master/readme-examples/nyquistplot-example.png)

#### Dependencies

```
library(pracma)
library(ggplot2)
```
### Sigma plot

Take the matrices A, B, C, and D of a linear MIMO state space representation and generate a sigma plot, i.e. a plot of the maximum and minimum singular values of the transfer function matrix evaluated over the frequency range.

#### Usage

```
sigmaplot <- function(fr, A = NULL, B = NULL, C = NULL, D = NULL)
```

#### Example

In order to plot the state space represented by

```
                   / 0  1 |   1   0 \     
       / A | B \   |-1 -1 |   0   1 |
G(s) = |---|---| = |------|---------|
       \ C | D /   | 10 0 | 0.1   0 |
                   \  0 1 |   0 0.1 /
```

from 10^-2 to 10^2 the function inputs take the form

```
sigmaplot(fr = c(-2, 2), 
		   A = matrix(c(0, 1, -1, -1), nrow = 2),
		   B = diag(2), 
		   C = matrix(c(10, 0, 0, 1), nrow = 2), 
		   D = 0.1*diag(2))
```

and result in

![](https://github.com/pascalschulthess/freqr/blob/master/readme-examples/sigmaplot-example.png)

#### Dependencies

```
library(reshape2)
library(ggplot2)
```

## Contact
- [p.schulthess@lacdr.leidenuniv.nl](mailto:p.schulthess@lacdr.leidenuniv.nl)
- [pascalschulthess.de](pascalschulthess.de)
