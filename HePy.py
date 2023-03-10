# -*- coding: utf-8 -*-
"""
HePy 1.0.0

Modelling helium cooling ages in apatite and zircon.

Created on Sun Aug 8 2021

@author: pete33geo
https://github.com/pete33geo/HePy.py
"""
import sys
import numpy as np
from scipy import integrate
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.optimize import newton


class Model_He():
    
    def __init__(self,mineral,U,Th,r):
        """Model helium cooling ages in apatite and zircon, following the
        approach of Ketcham (2005). The Model_He() object holds grain-specific
        data, called when calculating model cooling ages
        
        Parameters
        ----------
        mineral : str
            'apatite' or 'zircon'
        U : int
            Total Uranium [ppm]
        Th : int
            Thorium [ppm]
        r : int
            Equivalent spherical radius [microns]
            
        Methods
        -------
        ejection :
            Calculates alpha ejection correction factor for each parent isotope
        solve : 
            Calculates a helium cooling age, given a time-temperature path
        test : 
            Tests apatite and zircon RDAAM functions
        """
        # Mineral-specific stopping distances from Ketcham et al. (2011)
        if mineral =='apatite':
            self.s238 = 18.81
            self.s235 = 21.80
            self.s232 = 22.25
        elif mineral =='zircon':
            self.s238 = 15.55
            self.s235 = 18.05
            self.s232 = 18.43
        else:
            print('\n Error: Mineral "apatite" or "zircon" must be selected \n')
            sys.exit()
        
        self.mineral = mineral
        
        # Decay constants (1/yr)
        self.lmbd238 = 1.55125e-10
        self.lmbd232 = 4.9475e-11
        self.lmbd235 = 9.8485e-10
        
        # Create model domain for alpha ejection and diffusion calculations
        self.r = r
        self.nodes = 513
        self.rstep = r/self.nodes
        self.steps = np.arange(1,self.nodes + 1)
        self.x = (self.steps - 0.5) * self.rstep
        
        # Convert U and Th to atoms/g
        self.U238 = ((U / 1e6) * 0.992745/238.02891) * 6.02214179e23
        self.U235 = ((U / 1e6) * 0.007204/238.02891) * 6.02214179e23
        self.Th232 = ((Th / 1e6) / 232.03805) * 6.02214179e23
        
        self.eject = self.ejection()

    
    def ejection(self):
        """Calculates alpha ejection correction factors for each parent isotope.
        Called internally by main .solve() function, or may be called to get
        alpha ejection and correction factors for applying to a real dataset.

        Returns
        -------
        eject : dict
            Alpha ejection (ej) and correction factor (ft).
        """        
        vol = self.x**3
        vol[1:] = np.diff(vol)
        
        def _alpha_eject(S,atom):
            """Calculate alpha ejection and correction factor(ft) for each
            isotope. Called for both AHe and ZHe data.
            """
            aej = atom * (0.5 + (((self.x**2 + self.r**2 - S**2) / (2 * self.x)) - self.x)/ (2 * S))
            aej[self.x < (self.r - S)] = atom
            
            d = np.sum((vol * aej))
            nd = atom * self.r**3
            ft = (d / nd)

            return aej,ft
        
        ej238,ft238 = _alpha_eject(self.s238,self.U238)
        ej235,ft235 = _alpha_eject(self.s235,self.U235)
        ej232,ft232 = _alpha_eject(self.s232,self.Th232)   
        
        self.eject = {'ej238' :ej238,
                      'ej235' :ej235,
                      'ej232' :ej232,
                      'ft238' :ft238,
                      'ft235' :ft235,
                      'ft232' :ft232 }

        return self.eject
    

    def _damage(self):
        """Calculates damage and annealing and returns diffusivity for each
           timestep.
        """
        if self.mineral =='apatite': # Flowers et al. (2009)
            C0 = 0.39528
            C1 = 0.01073
            C2 = -65.12969
            C3 = -7.91715
            alpha = 0.04672
            rmr0 = 0.83
            kappa = 1.04 - rmr0
        elif self.mineral =='zircon': # Guenthner et al. (2013)
            C0=6.24534
            C1=-0.11977
            C2=-314.937
            C3=-14.2868
            alpha=-0.05721

        dt = (self.t1 - self.t2) * 365.25 * 24 * 60 * 60
        dam = np.zeros((len(dt),len(dt)))

        def _d_anneal(dt,T_mean):
            """Calculate change in track length at each timestep"""
            d = ((C0 + C1 * ((np.log(dt) - C2) /
                             (np.log(1/T_mean) - C3)))**(1/alpha) + 1)**-1
            return d
    
        def _t_eqv(T_mean,d):
            """Calculate equivalent time"""
            t_eqv = (np.exp(C2 + (np.log(1 / T_mean) - C3) *
                            (((1 / d) - 1)**alpha - C0) / C1))
            return t_eqv

        #calculate track length for newly generated tracks
        new_tracks = np.diag_indices_from(dam)
        dam[new_tracks] = _d_anneal(dt,self.T_mean)
        
        #calculate annealing using 'equivalent time'
        for i in range(1,len(self.t)-1):
            teqv = _t_eqv(self.T_mean[i],dam[i-1,:i])
            teqv = teqv + dt[i]
            dam[i,:i] = _d_anneal(teqv,self.T_mean[i])
        
        #volume conversion
        if self.mineral =='apatite':
            dam[(dam>=rmr0)] = ((dam[(dam>=rmr0)] - rmr0) / (1 - rmr0))**kappa
            dam[(dam>=0.765)] = 1.6 * dam[(dam>=0.765)] - 0.6

            df = ((dam!=0.0) & (dam<0.765))
            dam[df] = 9.205 * dam[df] * dam[df] - 9.157 * dam[df] + 2.269
            
            self.damage = dam
        
        elif self.mineral =='zircon':
            dam = 1.25 * (dam - 0.2)
            dam[dam<(0.36 / 1.25 + 0.2)] = 0
            
            self.damage = dam


    def _rdaam_diffusivity(self):
        """Calculate diffusivity based on the radiation damage calculated in 
        the damage function. Following (Flowers et al. 2009) for apatite RDAAM,
        and Guenthner et al. (2013) for zircon RDAAM.
        """
        self._damage()
        
        R = 0.008314472
        
        product = lambda lmbd : np.exp(lmbd * self.t1) - np.exp(lmbd * self.t2)
        
        if self.mineral =='apatite':
            rhov = (8/8 * self.U238 * 3.19 * product(self.lmbd238)
                    + 7/8 * self.U235 * 3.19 * product(self.lmbd235)
                    + 6/8 * self.Th232 * 3.19 * product(self.lmbd232)
                    )
            
            # Parameters: Flowers et al. (2009)
            lambdaf = 8.46e-17
            lambdaD = 1.55125e-10
            eta = 0.91
            L = 8.15e-4
            
            anneal_d = self.damage * (lambdaf / lambdaD) * rhov * eta * L
            anneal_d = np.sum(anneal_d, axis=0)
            
            omega = 1e-22
            psi = 1e-13
            Do = 0.6071
            Ea = 122.3
            Et = 34
            
            trap_diff = psi * anneal_d + omega * anneal_d**3
            
            diffusivities =  1e8 * ((Do * np.exp(-Ea / (R * self.T_mean)))
                             / (trap_diff * np.exp(Et / (R * self.T_mean)) + 1))
            
            self.diffusivities = diffusivities
            
        if self.mineral =='zircon':
            alphai = (8 * self.U238 * product(self.lmbd238)
                      + 7 * self.U235 * product(self.lmbd235)
                      + 6 * self.Th232 * product(self.lmbd232)
                      )
                                  
            anneal_d = alphai * np.flip(self.damage)
            anneal_d = np.sum(anneal_d, axis=0)
            
            # Parameters: Guenthner et al. (2013)
            Ba = 5.48E-19
            SV = 1.669
            D0l = 193188
            El = 165
            D0N17 = 0.0034
            EaN17 = 71
            Lint_lattice = 45920
                    
            a = self.r * 1e-4
            fa = 1 - np.exp(-Ba * anneal_d)
            DI = 1 - np.exp(-Ba * anneal_d * 3)
            Lint = (4.2 / (fa * SV) - 2.5)
            Tau = (Lint_lattice / Lint)**2
            DTaua2 = (1 / Tau) * D0l * np.exp(-El / (R * self.T_mean)) / (a * (1 - DI))**2
            DN17a2 = D0N17 * np.exp(-EaN17 / (R * self.T_mean)) / (a * DI)**2
            
            diffusivities = self.r**2 * (DI / DN17a2 + (1 - DI) / DTaua2)**-1
            
            self.diffusivities = diffusivities

    
    def _thermal_diffusion(self):
        """Diffusivity calculated as a function of temperature, with no rdaam. 
        Zircon parameters from Reiners et al. (2004); Apatite parameters from
        Farley et al. (2000).
        """
        if self.mineral == 'apatite': # Farley et al. (2000)
            Do = 50
            Ea = 137.522
        if self.mineral == 'zircon': # Reiners et al. (2004)
            Do = 0.46
            Ea = 169.0336
        
        R = 0.00831447
        T = self.T + 273.15
            
        d = (Do * np.exp(-Ea / (R * T))) * 1e8
        d = (d[:-1] + d[1:]) / 2 #average for dt
        
        self.diffusivities = d


    def _diffusion_solve(self):
        """Calculates a model cooling age given grain data, diffusivity, and a 
        time-temperature path.
        
        Returns
        -------
        model_age : float
            Model age [Myrs]
        """
        product = lambda lmbd : np.exp(lmbd * self.t1) - np.exp(lmbd * self.t2)

        #production term A in atoms/g
        A = (8 * self.eject['ej238'][:,None] * product(self.lmbd238) 
             + 7 * self.eject['ej235'][:,None] * product(self.lmbd235)
             + 6 * self.eject['ej232'][:,None] * product(self.lmbd232))
                
        k = self.diffusivities

        def _1D_solver():
            """Solve the 1D diffusion eqn for a sphere, using crank-nicolson"""
            dt = (self.t1 - self.t2) * 365.25 * 24 * 60 * 60
            dr = np.full((self.nodes,1),self.rstep)
            dr[0] = dr[0] / 2
            
            beta = np.zeros((np.shape(A)))
            beta = beta + ((2 * dr**2) / (k * dt))
            
            ad0 = (-2 - beta)
            bd0 = (2 - beta)
            
            ad0[:,0] = - 3 - beta[:,0]
            bd0[:,0] = 3 - beta[:,0]

            kd1 = np.ones((len(beta)-1))
            
            u = np.full((self.nodes),100.0,dtype=np.float64)
            
            f = (A * beta * self.x[:,None])
            
            a = sparse.diags([kd1,ad0[:,0],kd1],[-1,0,1],format="csr")
            b = sparse.diags([-kd1,bd0[:,0],-kd1],[-1,0,1],format="csr")
            
            for i in range(len(self.t1)):      
                a.setdiag(ad0[:,i])
                b.setdiag(bd0[:,i])
                u = spsolve(a, b.dot(u) - f[:,i])

            return u

        u = _1D_solver()

        integral = (u / self.x) * 4 * np.pi * self.x**2
        
        totalHe = integrate.romb(integral,dx=self.rstep)
        totalHe = totalHe / (4 * np.pi / 3)
        
        total238 = (8 * self.U238 * self.r**3) * self.eject['ft238']
        total235 = (7 * self.U235 * self.r**3) * self.eject['ft235']
        total232 = (6 * self.Th232 * self.r**3) * self.eject['ft232']
        
        He = total238 + total235 + total232 + totalHe
        
        m_He = lambda age: He - (total238 * np.exp(self.lmbd238 * age)
                                 + total235 * np.exp(self.lmbd235 * age)
                                 + total232 * np.exp(self.lmbd232 * age))
                
        model_age = newton(m_He,self.t[0],rtol=100) / 1e6
        
        return model_age

    
    def solve(self,t,T,k='rdaam'):
        """Main function to calculate He age for an inputted tT path.
        Diffusivity may use RDAAM (Apatite = Flowers et al. [2009];
        Zircon = Guenthner et al. [2013]), or simple diffusion
        (Apatite = Farley et al. [2000]; Zircon = Reiners et al. [2004]
        with no damage or annealing. Inherits alpha ejection correction.

        Inputs
        ------
        t : np.array
            Time [Myrs], resampled so max(dt) < 0.5
        T : np.array
            Temperature [C]
        k : str (optional)
            Diffusion model 'rdaam' (default) or 'simple'

        Returns
        -------
        He_age : float
            Model He age [Ma]
        """
        if (len(t)!=len(T)):
            print('\n Error: t and T arrays must be equal in size \n')
            sys.exit()
        if np.all(np.diff(t) > 0):
            print('\n Error: t values must decrease monotonically \n')
            sys.exit()
        
        self.T = T
        self.t = t
        self.t1 = t[:-1] * 1e6
        self.t2 = t[1:] * 1e6
        self.T_mean = (T[1:] + 273.15 + T[:-1] + 273.15) / 2

        if k == 'rdaam':
            self._rdaam_diffusivity()
        elif k == 'simple':
            self._thermal_diffusion()
        else:
            print('\n Error: k_model rdaam or diff must be selected \n')
            sys.exit()
        
        He_age = self._diffusion_solve()

        return He_age
    
    
    def test(self):
        """Test apatite and zircon RDAAM functions"""
        t_test = np.array([100,66,22,11,0])
        T_test = np.array([240,200,140,60,0])
        
        # Test zircon
        zr_test = Model_He(mineral='zircon',U=10,Th=40,r=100)
        zr_result = zr_test.solve(t_test,T_test,k='rdaam')
        zr_expected = 20.1202
        
        if np.allclose(zr_result,zr_expected,rtol=1e-05):
            print("Zircon RDAAM solve passed")
        else:
            print("Error: Zircon RDAAM solve failed")
        
        # Test apatite
        ap_test = Model_He(mineral='apatite',U=10,Th=40,r=100)
        ap_result = ap_test.solve(t_test,T_test,k='rdaam')
        ap_expected = 10.5859
        
        if np.allclose(ap_result,ap_expected,rtol=1e-05):
            print("Apatite RDAAM solve passed")
        else:
            print("Error: Apatite RDAAM solve failed")
