# HZ-reconstruction

Fortran code for constrain cosmology model using the 6 time delay measurements from H0LiCOW

To be used with the CosmoMC package ( https://cosmologist.info/cosmomc/ )

Note:
1. We take a snapshot of the original time delay distance chains from https://shsuyu.github.io/H0LiCOW/site/h0licow_data.html or https://github.com/shsuyu/H0LiCOW-public/tree/master/h0licow_distance_chains and they are stored under the direcory “original MCMC chains”.

2. In order to also use these time delay measurements in the CosmoMC, we have converted these chains into likelihood ( in the likelihood file).

3. Rsld.f90 for B1608 and the Rtdl.f90 for the other five time delay measurements. 


For more information, please refer to https://arxiv.org/abs/2001.08713
