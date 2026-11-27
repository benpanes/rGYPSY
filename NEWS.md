# rGYPSY 1.0.0

First release.

- `GYPSY_Sim()` projects stands year by year with GYPSY 2009. It follows the
  GYPSY SAS program (GYPSY_run.sas) and, where the GYPSY 2009 DLL differs
  from the SAS, the DLL.
- The model's components are exported: `GYPSY_SI()`, `GYPSY_SI_convert()`,
  `GYPSY_Age()`, `GYPSY_Ht()`, `GYPSY_PSI()`, `GYPSY_SDF()`, `GYPSY_N()`,
  `GYPSY_BAinc()` and `GYPSY_Vol()`.
- `GYPSY_Sim_TI()` and `GYPSY_SI_TI()` adjust site index for improved white
  spruce and lodgepole pine seedlots (Luo and Thomas 2021).
- SI, age and density types can be given per row.
- The model is a header-only C++ core with an Rcpp wrapper, shared with
  webGYPSY.
- Warnings for inputs outside GYPSY's range: total ages below the minimum,
  aspen younger than its years to breast height, low densities, and basal
  area given for species below breast height.
