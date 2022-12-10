import numpy as np
from scipy.integrate import quad
from scipy.special import spherical_jn

def no_wiggles_EH(cosmo_camb,z,k):
	cosmo=cosmo_camb
	T_cmb=cosmo.TCMB/2.7
	h=cosmo.h
	Om_0=cosmo.omch2+cosmo.ombh2
	Ob_0=cosmo.ombh2

	k_eq=7.46*Om_0/(T_cmb*100)**2
	sh=cosmo.h * 44.5 * np.log(9.83/Om_0)/np.sqrt(1 + 10 *Ob_0** 0.75)
	k=p.array(k)
	k=k[k>0]
	ks=k*sh
	q=k/(13.41*k_eq)
	alpha=1-0.328*np.log(431*Om_0)*(Ob_0/Om_0)+.38*np.log(22.3*Om_0)*(Ob_0/Om_0)**2
	gamma_eff=Om_0/h*(alpha+((1-alpha)/(1+(.43*k*sh)**4)))
	q_eff=Om_0*q/gamma_eff
	L0 = np.log(2*np.e + 1.8 * q_eff)
	C0 = 14.2 + 731.0 / (1 + 62.5 * q_eff)
	T= L0 / (L0 + C0 * q_eff**2)
	return T * cosmo.scale_independent_growth_factor(z)

def Sigma2_integrand(k,cosmo,z,bao_scale):
	Pnw=no_wiggles_EH(cosmo,z,k)
	integrand=Pnw*(1-spherical_jn(0,k*bao_scale)+2*spherical_jn(2,k*bao_scale))
	return integrand
	
def delta_Sigma2_integrand(k,cosmo,z,bao_scale):
	Pnw=no_wiggles_EH(cosmo,z,k)
	integrand=Pnw*(spherical_jn(2,k*bao_scale))
	return integrand
	
def sigma_tot2(mu,cosmo,q,bao_scale,k_s,z):
	'''
	mu: 
	k: wave-number (1/Mpc)
	bao_scale: bao scale in Mpc
	k_s: separation scale controling the modes (A.7 https://arxiv.org/pdf/1909.05277.pdf)
	Pnw: No wiggles power spectrum function
	'''
	f=((cosmo.omch2+cosmo.ombh2)/cosmo.h**2)**0.55 #flat cosmology
	
	
	Sigma2= quad(Sigma2_integrand,0,k_s,args=(bao_scale,cosmo,z))/(6*np.pi**2)
	delta_Sigma2=quad(delta_Sigma2_integrand,0,k_s,args=(bao_scale,cosmo,z))/(2*np.pi**2)
	sigma_par2=Sigma2*(1+f)**2
	return (mu**2)*sigma_par2+(1-mu)*Sigma2+(f*mu**2)*((mu**2) - 1)*delta_Sigma2
	
	
