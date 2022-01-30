# Model Selection for Scalar-on-Function Regression with Applications to Near-Infrared Spectroscopy

This is the final project for the **Research Module in Econometrics and Statistics** (Winter Semester 2021/2022) of [Jonghun Baek](https://github.com/Hun-baek), [Jakob R. Juergens](https://jakobjuergens.com/), and [Jonathan Willnow](https://github.com/JonathanWillnow).
This course is part of the Master of Science Economics at the University of Bonn and was supervised by Prof. Dr. Dominik Liebl and Dr. Christopher Walsh.

This repository contains the [Paper](Paper/main.pdf), the [slides for the presentation](Presentation/slides.pdf), and the [code](code) used for the simulations.

## Overview
In this project we studied model selection using cross validation in a scalar-on-function regression setting. The theoretical part of this paper introduces the necessary theory to conduct scalar on function regression using both basis expansions and functional principal components. The practical part, split into a Monte-Carlo simulation and an applciation to Near-Infrared Spectroscopy Data, is inspired by a paper called [**Functional Principal Component Regression and Functional Partial Least Squares**](https://www.tandfonline.com/doi/abs/10.1198/016214507000000527) by Philip T. Reiss and R. Todd Ogden.

In the simulation part we study the three different bases: b-splines, fourier and monomial and compare their performance for different numbers of basis functions in a cross-validation procedure focussed on prediction accuracy for a known coefficient function. Therefore, we do estimation using both a simple basis expansion of the observations and the coefficient function with the chosen basis, but also an estimation procedure based on the scores of functional principle components calculated based on the chosen basis.

The results of these simulations then inform the interpretation of the application where we try to predict real octane values from near-infrared spectroscopy curves from gasoline samples.
