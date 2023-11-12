# PeriodicWW_1D

Calculation of periodic non-linear water waves over finite depth using 
the Hamiltonian Coupled-Mode System (HCMS)

The function mainPropagation calculates the evolution of a 2nd order Stokes wave

References:

1) Papoutsellis (2016) [*Nonlinear water waves over varying bathymetry: Theoretical and numerical study using variational methods*](https://dspace.lib.ntua.gr/xmlui/handle/123456789/44741?locale-attribute=en) (Thesis NTUA 2016)

2) Papoutsellis & Athanassoulis (2017) [*A new efficient Hamiltonian approach to the nonlinear water-wave problem over arbitrary bathymetry*](https://arxiv.org/abs/1704.03276)

3) Athanassoulis & Papoutsellis (2017) [*Exact semi-separation of variables in waveguides with non-planar boundaries*](https://royalsocietypublishing.org/doi/10.1098/rspa.2017.0017), Proc. R. Soc. A **473**: 20170017.http://dx.doi.org/10.1098/rspa.2017.0017
   (Analytical expressions of the derivatives of the local wave-numbers)

4) Papoutsellis, Charalampopoulos, Athanassoulis (2018) [*Implementation of a fully nonlinear Hamiltonian Coupled-Mode Theory, and application to solitary wave problems over bathymetry*](https://www.sciencedirect.com/science/article/pii/S0997754617305617) Eur. J. Mech. B Fluids **72** 199–224
   (Analytical expressions of the matrix coefficients)

5) Papathanassiou, Papoutsellis, Athanassoulis (2019) [*Semi-explicit solutions to the water-wave dispersion relation and their role in the non-linear Hamiltonian coupled-mode theory*](https://link.springer.com/article/10.1007/s10665-018-09983-1) J. Eng. Math. **114**, 87–114
   (Semi-explicit calculation of local wave numbers and application to strongly non-linear periodic travelling waves)



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

