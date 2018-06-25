# Dual Approximate Dynamic Programming (DADP) for Stochastic Unit Commitment
Authors: Jagdish Ramakrishnan, James Luedtke

This repo contains the code and data for the paper: http://www.optimization-online.org/DB_HTML/2018/06/6672.html . The code is written in Matlab and optimization problems are solved using CPLEX 12.6 through the GAMS modeling language.

## Generator data instances

The data is stored in the .gdx file format for ease loading into GAMS. The files can be found in the /data generator directories. The instances have the form

uc{1}\_{2}md\_{3}sig.gdx,

where {1} is the number of generators, {2} is a scaling factor for the mean of the demand scenarios (\mu * 10 in the paper), i.e., {2} * TotCap/10 is the mean, and {3} is the percentage variation in the mean demand (\sigma in the paper).

## Running the code

### Example run

This would run all algorithms (Adelman \lambda_t approach, DADP approach, and perfect information approach) for the data instance uc15_4md_15sig.gdx. Enter into Matlab console, make sure all data files are in directory, and run this command:

```bash
>>> lb_main('uc15_4md_15sig.gdx',250,50/(4 * 0.1 * 4660),(0.99).^(0:249))
```

The above runs stochastic subgradient method with maximum iterations of 250, a constant gradient multiplier of 50/(4 * 0.1 * 4660), and decaying learning rate of 0.99. Note here TotCap = 4660.

### Short description of the main files

lb_main.m : main file which runs all three algorithms to generate bounds (i.e., Adelman \lambda_t approach, DADP \lambda_t(v_t) approach, and perfect information bound). Need to provide instance filename, max iterations / \rho constant for stochastic subgradient method, and step size vector.

max_Lg_r.m : algo for state independent Adelman \lambda_t approach

max_Lg_r_d.m : algo for DADP approach

run_per_info_mip.m : runs perfect information approach MIP

## Contact

If you have questions, feel free to reach out to me. My contact details are available at https://jramak.github.io/
