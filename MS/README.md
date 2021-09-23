# Critical *p* for Mølmer-Sørensen random unitaries

In order to find *p<sub>c</sub>* for MS systems, we need to run PastaQ simulations of these systems over many values of *N*, *&Theta;*, and *p*: these are all performed on Graham, which calls `ms_rand.sh` as an array job.  `ms_rand.sh` runs `ms_rand.jl`, which takes a randomly-sampled but pre-set set of `(N, Θ, p)` tuples.  This is important because I wanted to minimize the variation in job runtime for the different array jobs so that Graham can perform optimally, and expected runtimes increase with all three of the parameters.  This is still not quite ideal, as sometimes the jobs are not entirely finished because of the stochastic nature of the runtime, but I think this is pretty close to ideal: the runtimes are distributed roughly like a decaying exponential, and some will always slip through.  I began to implement a more sophisticated system by which the parameters could be assigned more dynamically, with each running job checking some communal file to see whatever parameters need to be computed at runtime, but this proved a little difficult because Graham doesn't like it when several jobs try to open the same file simulataneously, and I figured that the current solution was good enough.

After all of these data are computed by `ms_rand`, it needs to be analyzed to find the *p<sub>c</sub>*.  This is performed by `avg.jl`.  It does the same procedure as detailed in [Czischek 2021](https://arxiv.org/abs/2106.03769).  We expect our entropies to agree with 

<img src="https://user-images.githubusercontent.com/5233686/130818493-5d718932-86ea-47ff-81c4-73aa00242b6f.png" alt="drawing" width="400"/>

so for each potential value of *p<sub>c</sub>*, we perform a linear regression of the logarithms of each side of that equation.  We then just find the value of *p<sub>c</sub>* for each *&Theta;* that minimizes the sum-of-squares residual.  Here is a plot of the residuals that we have for each *&Theta;* as a function of *p<sub>c</sub>*.

![image](https://user-images.githubusercontent.com/5233686/130820506-97f97b41-3045-462b-ae48-9d3ae2369241.png)

Above is the data from late August—it is clearly not very good.  I had run about 800 jobs so far, but of the 2016-element set of parameter combinations, the minimum number of runs was 34, the maximum was 107, and the average was 98.  This is not enough to see clean data yet, so I sent Graham to run more jobs.


![image](https://user-images.githubusercontent.com/5233686/134468282-5fd2b987-76de-419a-aa8d-eb2372a76f61.png)

I set Graham to run this more times on August 23, and it is still running today, September 22.  Above is the current data.  The improvement is not great.  At this point, the min-ave-max for runs is 41-111-122.  The local minima are not really visible, and it's not clear that they'll necessarily emerge with more data, but it's really still too early to tell here.

# Code Overview

Here are the most relevant files in this directory and brief descriptions of their purpose:
- `ms_rand.jl` is the Julia file that actually runs these simulations.  It serializes its results into `msdata`.
- `ms_rand.sh` is the shell file that calls an array job with SLURM.  In order to run more jobs, the bounds that describe its array start and stop should be adjusted so that it won't overwrite any previous `msdata` files.
- `avg.jl` takes all of the current serial files in `msdata` and creates the serial files `Rs` and `denoms`.
- `Rs` is what is plotted above: it's the residuals of the fit of our data to the critical exponent expression.
- `denoms` is an array of how many runs each parameter tuple has had.  Its name comes from the fact that its values are the denominators when computing average entanglement.
