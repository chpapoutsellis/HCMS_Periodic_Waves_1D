# PeriodicWW_1D

Calculation of periodic non-linear water waves over finite depth using 
the Hamiltonian Coupled-Mode System (HCMS)

The function mainPropagation calculates the evolution of a steady travelling wave
or a perturbed wave 

Details on the numerical scheme can be found in:

Papoutsellis, thesis

Any comments and questions are welcome!


# Working Example
[comment]: <> Running the [ITCMS.m](ITCMS.m) function in the Matlab Command Window (type ITCMS and press Enter) produces the IT corresponding to the $M_2$ tidal constituent with constant stratification $N=0.0015$ (1/s), Coriolis frequency $f = 0.0001$ (1/s). The amplitude of the horizontal barotropic current at infinity is $U_0 = 0.04$ (m/s). The depth at infinity is $3000$ m and the ridge has a criticality $0.8$ and relative height $0.5$. The calculated energy conversion rate is $1577.26$ (W/m) per unit ridge length. We give below typical plots (baroclinic stream function and horizontal velocity) and a video produced by the code:
[comment]: <> ![alt text](https://github.com/ChPapoutsellis/InternalTidesCMSv1.0/blob/main/OUTPUT/psi.png?raw=true)
[comment]: <> ![alt text](https://github.com/ChPapoutsellis/InternalTidesCMSv1.0/blob/main/OUTPUT/u.png?raw=true)
[comment]: <> ![alt text](https://github.com/ChPapoutsellis/InternalTidesCMSv1.0/blob/main/OUTPUT/VIDEO.gif?raw=true)

[comment]: <> # Acknowledgements
[comment]: <> For the colormaps we use the function cmocean.m (also provided here) written by Chad A. Greene of the Institute for Geophysics at the 
[comment]: <> University of Texas at Austin (UTIG), June 2016, using colormaps created by Kristen
[comment]: <> Thyng of Texas A&M University, Department of Oceanography.

[comment]: <> Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True 
[comment]: <> colors of oceanography: Guidelines for effective and accurate colormap selection. 
[comment]: <> Oceanography 29(3):9-13, http://dx.doi.org/10.5670/oceanog.2016.66.

