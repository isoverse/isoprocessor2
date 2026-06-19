# Package options

These options are best set via `ip_options()` and queried via
`ip_get_option()`. However, the base functions
[`options()`](https://rdrr.io/r/base/options.html) and
[`getOption()`](https://rdrr.io/r/base/options.html) work as well but
require an `isoprocessor2.` prefix (the package name and a dot) for the
option name. Setting an option to a value of `NULL` means that the
default is used. `ip_get_options()` is available as an additional
convenience function to retrieve a subset of options with a regular
expression pattern.

## Usage

``` r
ip_options(...)

ip_get_options(pattern = NULL)

ip_get_option(x)
```

## Arguments

- ...:

  set package options, syntax identical to
  [`options()`](https://rdrr.io/r/base/options.html)

- pattern:

  to retrieve multiple options (as a list) with a shared pattern

- x:

  name of the specific option to retrieve

## Value

`ip_options()` and `ip_get_options()` return a named list of option
values; `ip_get_option()` returns the value of the single requested
option.

## Functions

- `ip_options()`: set/get option values

- `ip_get_options()`: get a subset of option values that fit a pattern

- `ip_get_option()`: retrieve the current value of one option (option
  must be defined for the package)

## Options for the isoprocessor2 package

- `debug`: turn on debug mode

- `auto_use_ansi`: whether to automatically enable correct rendering of
  stylized (ansi) output in HTML reports from notebooks that call
  [`library(isoprocessor2)`](https://isoprocessor2.isoverse.org/). Can
  be turned off by calling
  `isoprocessor2::ip_options(auto_use_ansi = FALSE)` **before** calling
  [`library(isoprocessor2)`](https://isoprocessor2.isoverse.org/).

## Examples

``` r
# All default options
ip_get_options()
#> $debug
#> [1] FALSE
#> 
#> $auto_use_ansi
#> [1] TRUE
#> 
```
