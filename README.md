# 3D, 2D and 1D physical benchmark catchment models

:warning: **This repository is still under construction.**

## About the repository

This GitHub repository includes implementation of physical benchmark catchment models formulated and used in three-part paper "On the development and analysis of coupled
surface-subsurface models of catchments" by P. Morawiecki and P.H. Trinh. Two main models are implemented:
* 3D/2D hillslope model - describes evolution of pressure head $h_g(x,y,z,t)$ and surface water height $h_s(x,z,t)$ in a 3D hillslope in a tilted V-shape catchment. Implemented using Finite Volume Method. The model is presented in "Part 2. A three-dimensional benchmark model and its properties" [[1]](#1).
* 1D hillslope model - describes evolution of total groundwater and surface water depth $H_g(x,t)$ in hillslopes characterised by a thin soil layer. It is implemented numerically using 1D PDE solver and a range of analytic and semi-analytic approximations are available. The model is presented in "Part 3. Analytical solutions and scaling laws" [[2]](#2).

Apart from implementation of above models, the repository also includes codes that can be used to recreate figures from two papers cited above. The visual appearance of the figures is different than in paper, since for the published version the output data were illustrated using LaTeX TikZ library.

## Repository content

When finalised the repository will include following directories:

* MODELS - implementation of 2D/3D and 1D model with additional functions
* RESULTS - directory for storing output datasets
* FIGURES - directory for storing output figures

The main directory includes the following examples:

| Filename | Short description | Figures recreated |
| ---      | ---       | ---       |
| 3D_example.m | compute and visualise 3D model solution | Fig 6 (Part 2) |
| 3D_to_2D_verification.m | compute error of 2D approximation of 3D model | Fig 7-10 (Part 2) |
| 2D_example.m | compute and visualise 2D model solution | Fig 10 (Part 2) |
| 2D_sensitivity_analysis.m | run sensitivity analysis for eight parameters of 2D model | Fig 11 (Part 2) |
| 2D_Maxwell_scenarios.m | compare 2D model to results from intercomparison study [[3]](#3) | Fig. 13-14 (Part 2) |
| 1D_examples.m | present short- and long-time behaviour of 1D model | Fig. 2-3 (Part 3) |
| 1D_parameter_variation.m | present impact of parameters on steady state and hydrograph for 1D model | Fig. 4 (Part 3) |
| 1D_steady_state.m | compare 1D model steady state to its asymptotic approximation | Fig. 5 (Part 3) |
| 1D_saturated_zone.m | compute growth of saturated zone in 1D model | Fig 6 (Part 3) |
| 1D_characteristic.m | show characteristic diagram for satruated zone of 1D model | Fig. 7 (Part 3) |
| 2D_1D_approximations.m | compare 2D and 1D model, and their analytic approximations  | Fig. 8 (Part 3) |
| 1D_sensitivity_analysis.m | run sensitivity analysis for seven parameters of 1D model | Fig. 9 (Part 3) |
| 1D_vertical_profile.m | present dependence of $h_g$, $\theta$ and $f$ on depth | Fig. 11 (Part 3) |
| 1D_steady_state_error.m | compute error of leading order approximation for initial steady state | Fig. 12 (Part 3) |
| 1D_saturated_zone_error.m | compute error leading order approximation for satruated zone size | Fig. 13 (Part 3) |

## Licence

## References

<a id="1">[1]</a> P. W. Morawiecki and P. H. Trinh. On the development and analysis of coupled surface-subsurface models of catchments. Part 2. A three-dimensional benchmark model and its properties. T.B.C., 2022.

<a id="2">[2]</a> P. W. Morawiecki and P. H. Trinh. On the development and analysis of coupled surface-subsurface models of catchments. Part 3. Analytical solutions and scaling laws. T.B.C., 2022.

<a id="3">[3]</a>  Maxwell, Reed M., et al. "Surface???subsurface model intercomparison: A first set of benchmark results to diagnose integrated hydrology and feedbacks." Water resources research 50.2 (2014): 1531-1549.
