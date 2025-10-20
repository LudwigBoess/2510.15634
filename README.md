# Simulating the LOcal Web (SLOW) - IV:
## γ-ray Emission in the Local Universe

In this repository you will find all scripts and dependencies to reproduce the figures presented in [Böss et. al. (2025)](https://arxiv.org/abs/2510.15634).

These scripts require Julia >= 1.7 to be installed. I reccomend using v1.12.

To initialize the dependencies run `julia build.jl`.

This will install all packages as well as set up the folder structure.

For the simulation snapshots send an email to `lboess@usm.lmu.de` and I will provide the relevant cutouts of the simulation domain to you.

You can get all plots by running from the main repository directory: `julia src/Fig01.jl`.

If you use an interactive session just ignore the warnings about redefinition of global variables.