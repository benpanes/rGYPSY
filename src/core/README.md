# GYPSY core (C++)

The whole GYPSY model in plain, header-only C++17. Nothing here depends on R
or Rcpp, so the same code can be compiled for other front ends, for example to
WebAssembly with Emscripten. The R package's bindings are in `../rgypsy.cpp`.

| Header                 | Contents                                                        |
|------------------------|-----------------------------------------------------------------|
| `gypsy_eq.h`           | Component equations, SAS missing-value helpers, core types      |
| `gypsy_si.h`           | SI and age solvers                                              |
| `gypsy_stand.h`        | Percent stocking index and stand density factor for one stand   |
| `gypsy_project.h`      | Yearly projection                                               |
| `gypsy_sim.h`          | Full run: input checks, SI and age, PSI and SDF, projection     |
| `gypsy_components.h`   | Vector versions of each component (used by the R functions)     |

## Conventions

- Missing values are `NaN`, in and out.
- Species are `"AW"`, `"SB"`, `"SW"` and `"PL"` (index 0-3 in that order).
- Invalid inputs throw `gypsy::Error`. Warnings are returned as strings
  (`gypsy::Messages`) for the caller to show.
- Results are `gypsy::Table`s: ordered column names with one
  `std::vector<double>` per column.

## Running a simulation

```cpp
#include "gypsy_sim.h"

gypsy::SimInput in;
in.stand_label = {"1", "2"};          // name of each stand code
in.stand     = {0, 1, 1};             // stand code per row (0 .. n_stands - 1)
in.stand_age = {35, 48, 48};
in.spp       = {"PL", "AW", "SW"};
in.n         = {2400, 1100, 450};
in.si        = {gypsy::NaN, gypsy::NaN, gypsy::NaN};
in.age       = {35, 48, 31};
in.ht        = {13.5, 19.8, 9.6};
in.ps        = {gypsy::NaN, 70, 25};
in.ba        = {18.2, 21.5, 4.1};

gypsy::SimOptions opt;                // defaults match GYPSY_Sim() in R
opt.simlen = 150;

gypsy::SimResult res = gypsy::Sim(in, opt);
// res.rows: one row per stand and year; res.init: one row per stand;
// res.warnings: messages
```

Output rows are in stand code order, so give codes in the order you want the
stands projected. The SAS carries one value from the last year of a stand
into the first year of the next, so the order can matter.

## Building outside R

Any C++17 compiler works, for example:

```bash
g++ -std=c++17 -O2 -I src/core my_program.cpp -o my_program
```

For WebAssembly, compile a small wrapper that exposes `gypsy::Sim` with
Emscripten (for example with `embind`):

```bash
em++ -std=c++17 -O2 -I src/core --bind wrapper.cpp -o gypsy.js
```
