import numpy as np
from scipy import special.legendre as legen
from scipy.integrate import quad

def P_ell(mu,ell,pk_damp,cosmo,z1,z2,k,Sigma2):
	frac=(2*ell+1)/2
	
	return quad(pk_damp*legen(mu),-1,1,args=(cosmo,mu,z1,z2,k,Sigma2))
