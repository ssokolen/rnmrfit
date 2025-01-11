# rnmrfit

** The following R package has been deprecated in favour of a standalone Rust implementation available under the master branch. **

This package implements NMR lineshape fitting using the real and imaginary components of the data in the frequency domain. The core of the algorithm is built around the NLOPT nonlinear optimization library with a number of helper scripts designed to facilitate working with NMR data.

More information can be found in the accompanying article: https://doi.org/10.1016/j.jmr.2018.11.004

## Recent updates

Simple lower and upper bound generation has been cleaned up with the introduction of `set_absolute_bounds` and `set_relative_bounds`. The default bounds implemented through`set_conservative_bounds` have also been re-thought and should hopefully lead to better fits out of the box. More details in the examples below.

Basic JCAMP data import has been merged into the master branch, but more testing is needed to account for various JCAMP variations. Import is currently limited to JCAMP files that store both real and imaginary data in the frequency domain. Bug reports would be appreciated if anyone encounters issues. 

## Installation

The `rnmrfit` package can be installed directly from GitHub using `devtools`:

```
#!R

library(devtools)
install_github('ssokolen/rnmrfit')
```

The package is still under active development and some functionality may change going forward. Further documentation will be added over time.  

## Tutorials

Tutorials have been prepared for typical use-cases when fitting 1D data. They are arranged in roughly increasing order of complexity. Most have been written with an R beginner in mind.

The code snippets assume that the package has been loaded `library(rnmrfit)`; however, any external packages are explicitly typed out in every snippet.

### Table of contents

1. Loading 1D data
2. Fitting one peak
3. Fitting multiple peaks
4. Extracting data from fit objects
5. One region, multiple spectra
6. Multiple regions, multiple spectra
7. Changing peak type
8. Including convolution in the fit


### Loading 1D data

As it stands, the `rnmrfit` package relies on previously processed data that it loads from the `pdata` folder of an NMR experiment. The following creates an `NMRData1D` object that is used to store NMR data.


```
#!R

nmrdata <- nmrdata_1d('/path/to/nmr/experiment/')
```

You can replace `'/path/to/nmr/experiment/'` with something like `'C:/Data/Experiments/10/'`. By default, `nmrdata_1d()` will load the smallest directory found in the `pdata` folder. You can override this with the `procs.number` argument.

```
#!R

nmrdata <- nmrdata_1d('/path/to/nmr/experiment/', procs.number = 999)
```

Outside of performing a lineshape fit, there is not much that you can do with the data alone -- you can filter it to a specified range of chemical shift values and plot it for a quick inspection. That said, you'll probably want to use something like Topspin for full featured data processing and visualization. The package comes with some example 13C data in the `nmr.chol` variable. 

```
#!R

# Filter nmr.chol to a range of 20-40 ppm and inspect the results
nmrdata <- filter_1d(nmr.chol, 20, 40)
plot(nmrdata)

# To see both the real and imaginary data
plot(nmrdata, components = 'r/i')
```

Note for R beginners: you can type `?filter_1d` to see the help file for the `filter_1d()` function. `plot()` is a little different as it can be used for many different objects, so you'll need to use `?plot.NMRData1D`.


### Fitting one peak

The following code will use the pre-loaded `nmr.chol` 13C experiment. See the above section for loading your own NMR data.

The fit process begins by the definition of a "scaffold" around the fit, which includes the number of peaks, peak type (Lorentz, Guass or Voigt), a description of the baseline, and a number of other options. This scaffold also serves as an initial guess for the fit process. The minimal amount of information required to generate a scaffold is a list of peaks and a set of data used to estimate initial guesses. The peak list can be defined using relatively standard NMR notation. A doublet at 120 ppm with a coupling constant of 30 Hz would be encoded as `'120 d 30'`. The `'d'` can be replaced by a number of common abbreviations including s, t, q, quint, sext, sept, oct, dd, dt, td, and tt. The latter 4 expect two coupling values, e.g., `'120 dd 30 30'`.

```
#!R

# Defining peak list, the list() command is optional if there is only 1 peak
peaks <- list('121.7 s')

# Although the fit can be performed on the whole data range, 
# it's far more efficient to filter the full data down to a range of interest
nmrdata <- filter_1d(nmr.chol, 120, 123)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata)

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
plot(fit)
```

The current default is to assume Lorentz peak shape, include phase correction, and fit a 3rd order baseline spline with 3 equidistant knots. These options are likely to work in a lot of the usual cases (although 3 knots may not be enough when considering a large range with a heavily fluctuating baseline). The package includes another version of the `nmr.chol` data called `nmr.chol.bad` where a rolling baseline has been added along with a 15 degree phase error (0th and 1st order). As can be seen in the code below, the default options have no problem with these errors.

```
#!R

# Defining peak list, the list() command is optional if there is only 1 peak
peaks <- list('121.7 s')

# Using the modified cholesterol data
nmrdata <- filter_1d(nmr.chol.bad, 120, 123)

# Examine the spectrum prior to fit
plot(nmrdata)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata)

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
plot(fit)
```

Although the default parameters work in this case, other cases may require some fine tuning. There are four optional parameters to consider: `include.phase` is a binary `TRUE`/`FALSE` toggle, `baseline.degree` is the polynomial order of each baseline spline, `n.knots` is the number of interior spline knots (the overall baseline is divided into `n.knots + 1` regions), and `include.difference` is another `TRUE`/`FALSE` toggle that controls whether separate baselines should be used for the real and imaginary data. Although both real and imaginary data should theoretically share the same baseline, difference could arise from a number of possible sources. The following is an example of how the default options can be overriden:

```
#!R

# Defining peak list, the list() command is optional if there is only 1 peak
peaks <- list('121.7 s')

# Using the modified cholesterol data
nmrdata <- filter_1d(nmr.chol.bad, 120, 123)

# Examine the spectrum prior to fit
plot(nmrdata)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata, include.phase = TRUE,
                           baseline.degree = 3, n.knots =3, 
                           include.difference = TRUE)

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
plot(fit)
```

If you check `?nmrfit_1d` you will see that by default, the fit function generates a set of conservative upper and lower bounds on the parameter estimates. It may necessary to to override these parameters by calling the boundary generation functions explicitly. There are two primary functions to set simple upper and lower boundaries: `set_absolute_bounds` and `set_relative_bounds`. "Absolute" bounds are applied the same to every single peak (such as a minimum and maximum peak width/height), whereas  "relative" bounds are applied in relation to the current value of the parameters (such as preventing the chemical shift position from changing by more than 1% of the fit data range starting from initial positions). Both absolute and relative bounds can set based on parameters that have been normalized to the observed data in the fit range (0-1 x-axis values and 0-1 y-axis values) or the original values in ppm/Hz/intensity. 

```
#!R

# Defining peak list, the list() command is optional if there is only 1 peak
peaks <- list('121.7 s')

# Using the modified cholesterol data
nmrdata <- filter_1d(nmr.chol.bad, 120, 123)

# Examine the spectrum prior to fit
plot(nmrdata)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata, include.phase = TRUE,
                           baseline.degree = 3, n.knots =3, 
                           include.difference = TRUE)

# Setting some manual constraints

# Limiting phase correction to pi/6 and setting peak width maximum (in Hz)
scaffold <- set_absolute_bounds(scaffold, phase = c(-pi/6, pi/6), 
			        width = c(0.3, 30), normalized = FALSE, 
				peak.units = 'hz')

# Preventing negative peak heights
scaffold <- set_absolute_bounds(scaffold, height = c(0, Inf), normalized = TRUE)

# Checking current bounds
print(lower_bounds(scaffold))
print(upper_bounds(scaffold))

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
plot(fit)
```
Note that tight constraints require good initial estimates. Keep in mind that small errors in initial position have a strong influence on initial peak height.


### Fitting multiple peaks

Fitting multiple peaks works exactly the same way as fitting one peak. However, it can be convenient to define a different name for each singlet or multiplet. Taking a look at a different region of the cholesterol data:

```
#!R

# Defining peak list (adding names is optional)
peaks <- list('C11'='21.116 s', 
              'C15'='24.339 s', 
              'C23'='23.892 s', 
              'C26'='22.645 s', 
              'C27'='22.915 s')

# Filtering to range of interest
nmrdata <- filter_1d(nmr.chol, 20.5, 26)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata)

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
plot(fit)
```

Taking a look at the resulting plot reveals another peak at approximately 22.78 ppm. To add it, it's necessary to generate a new scaffold.

```
#!R

# Defining peak list (adding names is optional)
peaks <- list('C11'='21.116 s', 
              'C15'='24.339 s', 
              'C23'='23.892 s', 
              'C26'='22.645 s', 
              'C27'='22.915 s',
              'unknown'='22.78 s')

# Filtering to range of interest
nmrdata <- filter_1d(nmr.chol, 20.5, 26)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata)

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
plot(fit)
```

By default, the `plot()` function combines all the lineshapes together. However, it may be usefult to plot them on their own. Although the default setting does not add the individual fit ids to the plot, this behaviour can be overridden with the `fit.legend` toggle.

```
#!R

plot(fit, sum.lineshapes = FALSE, legend.position = "right", fit.legend = TRUE)
```

See `?plot.NMRFit1D` for a couple more available options. Only a few display options have been added as it is expected that users would be more interested in making their own plots with the available data (as explained below).


### Extracting data from fit objects

Once a fit is generated, a number of parameters are of direct interest. This includes the calculated phase and lineshape parameters. Since baseline parameters are just terms of a spline function, they are more difficult to interpret directly. Use `phase()` and `peaks()` functions to extract the corresponding components of the overall fit. Using the fit from above:

```
#!R

# Defining peak list (adding names is optional)
peaks <- list('C11'='21.116 s', 
              'C15'='24.339 s', 
              'C23'='23.892 s', 
              'C26'='22.645 s', 
              'C27'='22.915 s',
              'uknown'='22.78 s')

# Filtering to range of interest
nmrdata <- filter_1d(nmr.chol, 20.5, 26)

# Generating scaffold using default parameters
scaffold <- nmrscaffold_1d(peaks, nmrdata)

# Fitting using default parameters
fit <- nmrfit_1d(scaffold)

# Examining fit
print(phase(fit))
print(peaks(fit))
```

As it stands, phase is limited to a single constant value for the entire region of the fit, but this will soon change. Note that the `width` column of the peaks may be in either Hz or ppm. It's possible to check and set the peak units using the `peak_units()` function.

```
#!R

# Check current peak units
print(peak_units(fit))

# Set to ppm and confirm change in width column
peak_units(fit) <- 'ppm'
print(peaks(fit))

# Set to hz and confirm change in width column
peak_units(fit) <- 'hz'
print(peaks(fit))
```

The peaks `data.frame` can be save to an excel readable form with the `write.csv()` command.

```
#!R

write.csv(peaks(fit), file = 'peak_data.csv')
```

Areas can be calculated for each multiplet using the `calc_area()` command, which outputs another `data.frame` that can also be saved.

```
#!R

areas <- calc_area(fit)
print(areas)
write.csv(areas, file = 'area_data.csv')
```

The exact lineshape data can be extracted using the `calc_lineshape()` function.

```
#!R

lineshapes <- calc_lineshape(fit)
print(head(lineshapes))
write.csv(lineshapes, file = 'lineshape_data.csv')
```

Note that the lineshape is complex (i.e. has real and imaginary components) and the peak data is stacked on top of each other. There are a number of ways of changing this, but the following example uses the `dplyr` and `tidyr` packages:

```
#!R

library(dplyr)
library(tidyr)

lineshapes.table <- lineshapes %>%
                      mutate(intensity = Re(intensity)) %>%
                      spread(id, intensity)
print(head(lineshapes.table))
write.csv(lineshapes.table, file = 'lineshape_table.csv')
```

The baseline can be extracted in a similar fashion using `calc_baseline` instead of `calc_lineshape`. By default, both functions use the same chemical shift values as the original data in the fit. However, a different range of chemical shifts could also be used. 

```
#!R

new.chemical.shift <- seq(22.5, 23.5, length.out = 200)
lineshapes <- calc_lineshape(fit, new.chemical.shift)
print(head(lineshapes))
```


### One region, multiple spectra

Although most operations are oriented around one region and one dataset, the package were designed with batch analysis in mind. Although it's possible to perform analysis of multiple spectra in many different ways, the `dplyr` package provides a convenient approach.

First, define a set of operations that will be performed on each spectrum as a `function()`. This function will take a filename (or path to an nmr folder) as an input and return the calculated peak areas as an output.

```
#!R

f_area <- function(filename) {

  # Inside, everything looks the same as before, but the operation are all based
  # on the "filename" variable, which will depend on function input
  nmrdata <- nmrdata_1d(filename)
  nmrdata <- filter_1d(nmrdata, 20.5, 26)
  
  # Defining peak list
  peaks <- list('C11'='21.116 s', 
                'C15'='24.339 s', 
                'C23'='23.892 s', 
                'C26'='22.645 s', 
                'C27'='22.915 s',
                'uknown'='22.78 s')

  # Generating scaffold using default parameters
  scaffold <- nmrscaffold_1d(peaks, nmrdata)

  # Fitting using default parameters
  fit <- nmrfit_1d(scaffold)
  
  # Calculating areas
  areas <- calc_area(fit)

  # The last variable in the function is returned as output
  areas
}
```

It's possible to run the function defined above on only one directory.

```
#!R

path <- 'inept_z/15'
f_area(path)
```

However, the idea is to define a vector of different paths and 'map' the function to each path. The following code snippet uses the `dplyr`, `purrr`, and `tidyr` packages to help with the mapping.

```
#!R

library(dplyr)
library(purrr)
library(tidyr)

paths <- file.path('inept_z', 15:19)
areas <- data_frame(filename = paths)
print(areas)

areas <- areas %>%
           mutate(area = map(filename, f_area)) %>%
           unnest(area)

print(areas)
```

As before, it's possible to spread the result out into a table and save it to file.

```
#!R

library(dplyr)
library(tidyr)

areas.table <- areas %>%
                 spread(id, area)
print(head(areas.table))
write.csv(areas.table, file = 'areas_table.csv')
```


### Multiple regions, multiple spectra

It's possible to expand the `f_area` function that was defined above to consider different regions (or other variables). To start, we can define two new inputs: a lower and upper boundary.

```
#!R

f_area <- function(filename, lower, upper, peaks) {

  nmrdata <- nmrdata_1d(filename)
  
  # The lower and upper boundaries now change with different inputs
  nmrdata <- filter_1d(nmrdata, lower, upper)

  # The peaks variable must also be provided in the input
  # (as it will change for different lower/upper bounds)
  scaffold <- nmrscaffold_1d(peaks, nmrdata)

  # Fitting using default parameters
  fit <- nmrfit_1d(scaffold)
  
  # Calculating areas
  areas <- calc_area(fit)

  # The last variable in the function is returned as output
  areas
}
```

Once the function is defined, it's necessary to make a `data.frame` that combines the files to open, the lower and upper bounds, and the peak definitions. It's easiest to define the files separately, and then combine all unique combinations of files and peak regions.

```
#!R

library(dplyr)

# First the paths
paths <- file.path('inept_z', 15:19)
paths <- data_frame(filename = paths)
print(paths)

# Then the regions and peak definitions
regions <- data_frame(lower = c(120, 20.5),
                      upper = c(123, 26),
                      peaks = list(list('C6'='121.7 s'),
                                   list('C11'='21.116 s', 
                                        'C15'='24.339 s', 
                                        'C23'='23.892 s', 
                                        'C26'='22.645 s', 
                                        'C27'='22.915 s',
                                        'uknown'='22.78 s')))
                                        
# Combine all unique combinations of files and regions
areas <- crossing(paths, regions)
print(areas)
```

The mapping looks a little different as we've gone from only one variable (filename) to four (filename, lower, upper, peaks).

```
#!R

library(dplyr)
library(purrr)
library(tidyr)

areas <- areas %>%
           mutate(area = pmap(list(filename, lower, upper, peaks), f_area)) %>%
           unnest(area)

print(areas)
```


### Changing peak type

The peak type of the fit can be changed with the `peak.type` argument of the `scaffold()` function. Currently, four peak types are supported: `'lorenz', 'gauss', 'voigt', 'pvoigt'`. Whereas 'pvoigt' uses a simple addition of 'lorenz' and 'gauss' lineshapes (assigning two different heights), 'voigt' is a convolution of the two lineshapes (where the height is determined by the 'lorenz' component and the widths are different).

Comparing 'voigt' and 'pvoigt' options:

```
#!R

peaks <- list('121.7 s')
nmrdata <- filter_1d(nmr.chol, 120, 123)

scaffold <- nmrscaffold_1d(peaks, nmrdata, peak.type = 'voigt')
fit <- nmrfit_1d(scaffold)

print(peaks(fit))
plot(fit)

scaffold <- nmrscaffold_1d(peaks, nmrdata, peak.type = 'pvoigt')
fit <- nmrfit_1d(scaffold)

print(peaks(fit))
plot(fit)
```

In both cases, the 'gauss' component is very small for this peak.

It's worth pointing out that 'voigt' and 'pvoigt' fits are less robust than the standard 'lorenz' one. As such, it may be helpful to get initial fit values using the 'lorenz' peak type and then switch the type to 'voigt'. This will conserve previously calculated parameters and use them as the initial values in the 'voigt' calculation.

```
#!R

peaks <- list('121.7 s')
nmrdata <- filter_1d(nmr.chol, 120, 123)

# First fit using the standard 'lorenz' type
scaffold <- nmrscaffold_1d(peaks, nmrdata)
fit <- nmrfit_1d(scaffold)

# Change the peak type of the fit and store as new scaffold
scaffold2 <- set_peak_type(fit, peak.type = 'voigt')
fit2 <- nmrfit_1d(scaffold2)

print(peaks(fit2))
plot(fit2)
```

Since the fit was good the first time, this step was not technically necessary for this peak.


### Including convolution in the fit

The default fit uses the spectra as it appears in topspin. However, it's also possible to fit a lineshape as a convolution of an ideal Lorenzian (or Voigt) lineshape and an apodization function. Since truncation can be seen as apodization with a rectangular step function, this fit can explicitly capture the effect of truncation (among other possible artefacts). To take the apodization into account, it's necessary to specify the convolution using the `set_convolution()` function. By default, `set_convolution()` will read the apodization and zero-fill parameters from the `procs` parameters (specified by `lb`, `gb`, `ssb`, `si`). Note, it's important to make sure that these are all set to zero in programs like Topspin if no apodization is actually applied. An `nmr.chol.trunc` variable has been provided with the package to demonstrate how this works.

```
#!R

# Looking at only a portion nmr.chol.trunc
nmrdata <- filter_1d(nmr.chol.trunc, 120, 123)

# View truncated data
plot(nmrdata)

# Load the convolution from internal procs data
nmrdata <- set_convolution(nmrdata)
```

After than, the fit can be performed as normal; however, there are still a few caveats. There could be problems if the initial peak position falls on a side lobe and including convolution considerably slows down the fit process, but it should be possible to speed this up in future updates.

```
#!R

peaks <- list('121.74 s')
scaffold <- nmrscaffold_1d(peaks, nmrdata)
fit <- nmrfit_1d(scaffold)

plot(fit)
```

It's possible to separate truncation from the ideal lineshape.

```
#!R

# The include.convolution option removes convolution from the fit
# and deconvolve.residual removes the lack of fit from the residual line.
plot(fit, include.convolution = FALSE, deconvolve.residual = TRUE)

# Note the much smaller width of the peak parameters
print(peaks(fit))
```

The exact same procedure can be applied to directly fit other types of apodization. Custom apodization functions can also be supplied, see `?set_convolution` for more details.
