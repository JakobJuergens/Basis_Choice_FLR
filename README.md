# Basis Choice for Scalar-on-Function Regression with Applications to Near-Infrared Spectroscopy

> **Archive notice.** This is an archived master's-level collaborative student project. It is preserved as historical material and as an example for students considering a research-module project. It appears in the [archive of earlier projects](https://jakobjuergens.com/previous-projects/) and is not part of Jakob R. Juergens's [current research portfolio](https://jakobjuergens.com/research/).

## Project context

- **Course:** Final project for the Research Module in Econometrics and Statistics
- **Institution:** University of Bonn
- **Completed:** 2022 (Winter Semester 2021/2022)
- **Authors:** [Jonghun Baek](https://github.com/Hun-baek), [Jakob R. Juergens](https://jakobjuergens.com/), and [Jonathan Willnow](https://github.com/JonathanWillnow)
- **Supervisor:** [Prof. Dr. Dominik Liebl](https://www.dliebl.com/)
- **Additional supervision credited in the archived paper:** Dr. Christopher Walsh

## Start here

- Read the [final paper](Paper/main.pdf) for the complete theory, simulation study, and empirical application.
- View the [presentation slides](Presentation/slides.pdf) for a shorter account of the project.
- Browse the principal [code directory](Code/), beginning with [Main_Simulation.R](Code/Main_Simulation.R) for the Monte Carlo study and [Main_Application.R](Code/Main_Application.R) for the gasoline application.
- See the [contextualized project archive](https://jakobjuergens.com/previous-projects/) for a retrospective summary alongside other earlier student work.

## Overview

The project studies scalar-on-function regression, where a scalar outcome is predicted from a functional covariate. It compares direct basis-expansion regression with functional principal component regression (FPCR). Fourier, B-spline, and monomial bases are used to represent the functional observations, and cross-validated mean squared prediction error guides the choice of basis size or truncation and the number of functional principal components.

The simulation study compares these approaches across smooth and more irregular coefficient functions and different noise settings. Its results illustrate bias-variance and numerical-stability tradeoffs rather than establishing a universal ranking of bases. The empirical application predicts octane values from 60 near-infrared spectroscopy curves. The paper reports FPCR with four functional principal components constructed from seven Fourier basis functions as the best specification considered, with direct B-spline basis regression performing similarly.

## Repository guide

- [`Paper/`](Paper/) contains the paper PDF, LaTeX source, tables, and bibliography files.
- [`Presentation/`](Presentation/) contains the slides PDF and its LaTeX source.
- [`Code/`](Code/) contains the R analysis. [`Main_Simulation.R`](Code/Main_Simulation.R) runs the simulation specifications, [`Main_Application.R`](Code/Main_Application.R) runs the empirical specifications, and [`Main_Stitching.R`](Code/Main_Stitching.R) aggregates partial simulation output. Supporting functions live alongside these scripts.
- [`Code/Results/`](Code/Results/) preserves the archived simulation and application outputs, while [`Graphics/`](Graphics/) contains figures used in the paper and presentation.
- [`Code/Archive/`](Code/Archive/) contains earlier exploratory scripts and notebooks; it is historical material rather than the natural starting point.

Together, the paper and code show how the project combines functional-regression theory with simulation design, cross-validated model selection, and an empirical spectroscopy application.

## Reproducibility note

The analysis is written primarily in R; the paper and slides use LaTeX. The main R scripts load `refund`, `MASS`, `tidyverse`, `fda`, `fpca`, and `caret`. They use paths relative to [`Code/`](Code/), so that directory is the apparent working directory. There is no single end-to-end execution script: [`Main_Simulation.R`](Code/Main_Simulation.R) and [`Main_Application.R`](Code/Main_Application.R) are the respective entry points, with [`Main_Stitching.R`](Code/Main_Stitching.R) used to combine partial simulation results.

The gasoline spectroscopy data are not stored as a standalone raw-data file in this repository. The scripts load the `gasoline` data supplied by the `refund` R package via `data(gasoline)`, so that package is also the external data source for a fresh run. The repository has an RStudio project file but no dependency lockfile, and package versions are not pinned. This archival documentation pass did not comprehensively modernize or rerun the analysis under current R and package versions; the presence of archived outputs should not be read as a current reproducibility guarantee.

## Retrospective

This repository remains a compact example of a collaborative research project combining methodological exposition, simulation design, model selection, and an empirical application.
