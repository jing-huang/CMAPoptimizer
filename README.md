# CMAPoptimizer
Monte Carlo simulated annealing to optimize the CMAP potentials in molecular mechanics force fields with reweighting calculation.

A python code, together with example input data, is provided to optimize the 2D CMAP potential in force fields. Monte Carlo simulated annealing is carried out using a generic target function with reweighting technique. This code was used to optimize the non-Gly, non-Pro CMAP potentials in the CHARMM36m protein force field to improve the conformational sampling of left helix. The release entitled "paper" is the exact code to derive C36m CMAP, while further improvement continues. The calculation of CMAP energy is fully coded in python, which might be useful for other force field development projects.

See J. Huang, S. Rauscher, G. Nawrocki, T. Ran, M. Feig, B. de Groot, H. Grubmueller and A. MacKerell Jr.
CHARMM36m: An Improved Force Field for Folded and Intrinsically Disordered Proteins
Nature Methods, 14, 71 (2017)

