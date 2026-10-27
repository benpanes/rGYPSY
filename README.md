# rGYPSY

An R and Rcpp implementation of the Alberta Growth and Yield Projection System
(GYPSY 2009). `GYPSY_Sim()` implements the core model logic through multiple subfunctions and projects stands year by year.

The model's components are exported and can be called individually:
`GYPSY_SI()`, `GYPSY_SI_convert()`, `GYPSY_Age()`, `GYPSY_Ht()`,
`GYPSY_PSI()`, `GYPSY_SDF()`, `GYPSY_N()`, `GYPSY_BAinc()` and `GYPSY_Vol()`.

## Requirements

- R (developed on R 4.4).
- R packages [Rcpp](https://cran.r-project.org/package=Rcpp) and
  [data.table](https://cran.r-project.org/package=data.table). These are
  installed automatically with rGYPSY.

## Installation

Install the prebuilt package from
[R-universe](https://benpanes.r-universe.dev). No compiler is needed on
Windows or macOS:

```r
install.packages("rGYPSY",
                 repos = c("https://benpanes.r-universe.dev",
                           "https://cloud.r-project.org"))
```

### From source

To install the latest version from GitHub, or on Linux, you need a C++
compiler:

- Windows: [Rtools](https://cran.r-project.org/bin/windows/Rtools/), the
  version that matches your R version.
- macOS: the Xcode command line tools (`xcode-select --install`).
- Linux: usually already installed (`g++`).

```r
install.packages("remotes")
remotes::install_github("benpanes/rGYPSY")
```

If you have a local clone of the repository, you can install from that
instead:

```r
remotes::install_local("path/to/rGYPSY")
```

## Input data

`GYPSY_Sim()` takes one row per stand and species. Pass each column as a
separate argument:

| Argument    | Description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| `stand`     | Stand identifier. Each species can appear only once per stand.              |
| `stand_age` | Stand age. Copied to the `standage` output column; not used in the model.  |
| `spp`       | Species group: `AW` (aspen), `SB` (black spruce), `SW` (white spruce) or `PL` (lodgepole pine). |
| `n`         | Density (stems/ha).                                                         |
| `si`        | Site index (m). Breast height age 50 by default (`sitype = "B"`), or total age 50 with `sitype = "T"`. |
| `age`       | Age (years). Total age by default (`agetype = "T"`), or breast height (`"B"`) or stump height (`"S"`) age. |
| `ht`        | Top height (m).                                                             |
| `ps`        | Percent stocking (optional, `NA` if unknown).                               |
| `ba`        | Basal area (m²/ha, optional, `NA` if unknown).                              |

Each row needs two out of `si`, `age` and `ht` in order to calculate the third:

- With `age` and `ht`, SI is computed from them. Total age must be at least
  7.5 years for AW and PL and 9.5 for SB and SW.
- A supplied `si` takes precedence over the computed SI.
- With `si` and `ht`, age is computed from them, and any supplied `age` is
  ignored.

Optional inputs:

- **Percent stocking (`ps`).** Converted to a percent stocking index, which
  adjusts the density and basal area increment models.
- **Basal area (`ba`).** Used to localize the predicted basal area to the
  observed value.

Set `psflag = 0` or `baflag = 0` to ignore these columns.

An example input file is included with the package:

```r
system.file("extdata", "example_input.csv", package = "rGYPSY")
```

```
stand,stand_age,spp,n,si,age,ht,ps,ba
1,35,PL,2400,,35,13.5,,18.2
2,48,AW,1100,,48,19.8,70,21.5
2,48,SW,450,,31,9.6,25,4.1
3,62,SB,2900,12.5,,13.2,,22.8
3,62,PL,600,16.0,,18.1,,6.3
...
```

The stands in the example file cover these cases:

- **Stand 1:** a pure pine stand.
- **Stand 2:** a mixedwood with percent stocking.
- **Stand 3:** SI and height given instead of age.
- **Stand 4:** all four species.
- **Stand 5:** no basal area.

Blank cells are read as missing.

## Workflow example

```r
library(rGYPSY)
library(data.table)

# 1. Read the input data (one row per stand and species)
inp <- fread(system.file("extdata", "example_input.csv", package = "rGYPSY"))

# 2. Project each stand for 150 years
out <- GYPSY_Sim(stand     = inp$stand,
                 stand_age = inp$stand_age,
                 spp       = inp$spp,
                 n         = inp$n,
                 si        = inp$si,
                 age       = inp$age,
                 ht        = inp$ht,
                 ps        = inp$ps,
                 ba        = inp$ba,
                 simlen    = 150)

# 3. Yearly output: one row per stand and year, species in columns
out[year %in% c(50, 100), .(stand, totage, ba_pl, tvol_all, mvol_all)]

# 4. Starting values per stand (SI, ages, PSI, SDF, initial density, ...)
attr(out, "stand_init")

# 5. Merchantable volume at the culmination of merchantable MAI
out[, .SD[which.max(mai_m_all)], by = stand][, .(stand, totage, mvol_all, mai_m_all)]
```

The projection starts at the fractional part of the leading species' total
age and runs for `simlen` years.

The output has these columns:

- `stand`, `year`, `totage` (leading species total age) and `standage`.
- For each species (`aw`, `sb`, `sw`, `pl`):
  - `totage_` and `bhage_`: total and breast height age.
  - `topht_`: top height.
  - `n_`: density.
  - `sc_`: species composition.
  - `ba_`: basal area (m2/ha).
  - `qmd_`: quadratic mean diameter (cm).
  - `tvol_` and `mvol_`: total and merchantable volume (m3).
  - `mn_`: merchantable stems.
  - `mai_t_` and `mai_m_`: total and merchantable MAI (m3/yr).
  - Conifer and stand totals of volume and MAI.

See `?GYPSY_Sim` for the full list.

### Merchantability standards

The default standard is a 0.3 m stump height, a 15 cm stump diameter outside
bark and a 10 cm top diameter inside bark. To change it, give one value for
all species or a value per species:

```r
out_1310 <- GYPSY_Sim(inp$stand, inp$stand_age, inp$spp, inp$n,
                      si = inp$si, age = inp$age, ht = inp$ht,
                      ps = inp$ps, ba = inp$ba, simlen = 150,
                      stump_dob = c(AW = 15, SB = 13, SW = 15, PL = 13),
                      top_dib   = 10)
```
