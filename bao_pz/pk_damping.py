import camb
impot numpy as np

class no_wiggles_EH:
	def __init__(self,cosmo_camb,z):
		self.cosmo=cosmo_camb
		self.z=z
		self.T_cmb=cosmo.TCMB/2.7
		self.h=cosmo.h
		self.Om_0=cosmo.omch2+cosmo.ombh2
		self.Ob_0=cosmo.ombh2

		self.k_eq=7.46*self.Om_0/(self.T_cmb*100)**2
		self.sh=cosmo.h * 44.5 * np.log(9.83/self.Om_0)/np.sqrt(1 + 10 *self.Ob_0** 0.75)
	def __call_(self,k):
		k=p.array(k)
		k=k[k>0]
		ks=k*self.sh
		q=k/(13.41*self.k_eq)
		alpha=1-0.328*np.log(431*self.Om_0)*(self.Ob_0/self.Om_0)+.38*np.log(22.3*self.Om_0)*(self.Ob_0/self.Om_0)**2
		gamma_eff=self.Om_0/self.h*(self.alpha+((1-self.alpha)/(1+(.43*k*self.sh)**4)))
		q_eff=self.Om_0*q/gamma_eff
		L0 = np.log(2*np.e + 1.8 * q_eff)
		C0 = 14.2 + 731.0 / (1 + 62.5 * q_eff)

		T= L0 / (L0 + C0 * q_eff**2)
        return T * self.cosmo.scale_independent_growth_factor(self.redshift)
