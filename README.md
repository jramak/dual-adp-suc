# Dual Approximate Dynamic Programming (DADP) for Stochastic Unit Commitment
Authors: Jagdish Ramakrishnan, James Luedtke

This repo contains the code and data for the paper: http://www.optimization-online.org/DB_HTML/2018/06/6672.html . The code is written in Matlab and optimization problems are solved using CPLEX 12.6 through the GAMS modeling language. 

## Generator data instances

The data is stored in the .gdx file format for ease loading into GAMS. The instances have the form

uc{1}\_{2}md\_{3}sig.gdx, 

where {1} is the number of generators, {2} is a scaling factor for the mean of the demand scenarios (\mu * 10 in the paper), i.e., {2} * TotCap/10 is the mean, and {3} is the percentage variation in the mean demand (\sigma in the paper).

## Contact

If you have questions, feel free to reach out to me. My contact details are available at https://jramak.github.io/
