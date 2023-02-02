# HePy_v0.1.py

import sys
import numpy as np
from scipy import integrate
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time as runningtime


class Model_He():
    
    def __init__(self,mineral,U,Th,r):
        """Calculate a model cooling age for input grain data, and a
        time-temperature path. Time temperature path is inputted to the
        Model_He.solve() function. Following the approch of Ketcham, 2005.
        
        Parameters
        ----------
        mineral:
            'apatite' or 'zircon'
        U:
            Total Uranium in ppm
        Th:
            Thorium in ppm
        r:
            Equivalent spherical radius in microns
        
        Examples
        --------
        """
        #mineral-specific stopping distances
        if mineral =='apatite':
            self.s238 = 18.81
            self.s235 = 21.80
            self.s232 = 22.25
        elif mineral =='zircon':
            self.s238 = 16.97
            self.s235 = 19.64
            self.s232 = 19.32
        else:
            print('\n Error: Mineral (apatite or zircon) must be selected \n')
            sys.exit()
        
        self.mineral = mineral
                
        self.lmbd238 = 1.55125e-10
        self.lmbd232 = 4.9475e-11
        self.lmbd235 = 9.8485e-10
        
        #create model domain for alpha eject and diffusion calculations
        self.r = r
        self.nodes = 513 #can be any number - reduce to speed up?
        self.rstep = r/self.nodes
        self.steps = np.arange(1,self.nodes+1)
        self.x = (self.steps - 0.5) * self.rstep
        
        #convert U and Th to atoms/g
        self.U238 = ((U / 1e6) * 0.992745/238.02891) * 6.02214179e23
        self.U235 = ((U / 1e6) * 0.007204/238.02891) * 6.02214179e23
        self.Th232 = ((Th / 1e6) / 232.03805) * 6.02214179e23
        
        self.eject = self.eject_factor()
    
    def eject_factor (self):
        """Returns alpha ejection correction factors for each parent isotope.
        Called internally by main modelling function.
            
        Returns
        -------
        Dictionary
        """    
        x = self.x
        r = self.r
        
        vol = x**3
        vol[1:] = np.diff(vol)
        
        def _alpha_eject(S,atom):
            """ Nested function to calculate alpha ejection and correction
            factor(ft) for each isotope. May be used for both AHe and ZHe data.
            """
            aej = (0.5 + (((x**2 + r**2 - S**2) / (2 * x)) - x)/ (2 * S))
            aej = atom * aej
            
            aej[x < (r - S)] = atom
            
            d = np.sum((vol * aej))
            nd = atom * r**3
            
            ft = d / nd
            
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
            
        Returns
        -------
        diffusivity np.1Darray [microns**2/s]
        """
        if self.mineral =='apatite':
            #damage and annealing
            C0 = 0.39528
            C1 = 0.01073
            C2 = -65.12969
            C3 = -7.91715
            alpha = 0.04672
            rmr0 = 0.83
            kappa = 1.04 - rmr0
        elif self.mineral =='zircon':
            C0=6.24534
            C1=-0.11977
            C2=-314.937
            C3=-14.2868
            alpha=-0.05721

        dt = (self.t1 - self.t2) * 365.25 * 24 * 60 * 60
        dam = np.zeros((len(dt),len(dt)))
        T_mean = self.T_mean

        def d_anneal(dt,T_mean):
            """Calculate change in track length at each timestep"""
            d = ((C0 + C1 * ((np.log(dt) - C2) /
                             (np.log(1/T_mean) - C3)))**(1/alpha) + 1)**-1
            return d
    
        def t_eqv(T_mean,d):
            """Calculate equivalent time"""
            t_eqv = (np.exp(C2 + (np.log(1 / T_mean) - C3) *
                            (((1 / d) - 1)**alpha - C0) / C1))
            return t_eqv

        #calculate track length for newly generated tracks
        new_tracks = np.diag_indices_from(dam)
        dam[new_tracks] = d_anneal(dt,T_mean)
        
        #calculate annealing using 'equivalent time' (see e.g., ketcham 2005)
        for i in range(1,len(t)-1):
            teqv = t_eqv(T_mean[i],dam[i-1,:i])
            teqv = teqv + dt[i]
            dam[i,:i] = d_anneal(teqv,T_mean[i])
        
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

    def rdaam_diffusivity(self):
        """Calculate diffusivity based on the damage calculated in the fission
        tracks fn.
        """
        self._damage()
        
        t1 = self.t1
        t2 = self.t2
        T_mean = self.T_mean
        
        U238 = self.U238 
        U235 = self.U235
        Th232 = self.Th232
        
        lmbd238 = self.lmbd238
        lmbd235 = self.lmbd235
        lmbd232 = self.lmbd232
        
        R = 0.008314472
        
        if self.mineral =='apatite':
            rhov = (((8 / 8) * (U238 * 3.19)
                     * (np.exp(lmbd238 * t1) - np.exp(lmbd238 * t2))) 
                    + ((7 / 8) * (U235 * 3.19)
                       * (np.exp(lmbd235 * t1) - np.exp(lmbd235 * t2)))
                    + ((6 / 8) * (Th232 * 3.19)
                       * (np.exp(lmbd232 * t1) - np.exp(lmbd232 * t2))))
            
            lambdaf = 8.46e-17
            lambdaD = 1.55125e-10
            etaq = 0.91
            L = 8.15e-4
            omega = 1e-22
            psi = 1e-13
            Do = 0.6071
            Ea = 122.3
            Et = 34

            anneal_d = self.damage * (lambdaf / lambdaD) * rhov * etaq * L
            anneal_d = np.sum(anneal_d, axis=0)
                    
            trap_diff = psi * anneal_d + omega * anneal_d**3
            
            diffusivities = ((Do * np.exp(-Ea / (R * T_mean)))
                             / (trap_diff * np.exp(Et / (R * T_mean)) + 1))
            
#            self.diffusivities = diffusivities * 1e4#cm2/s to microns2/s ####### CHECK UNITS!
        
        if self.mineral =='zircon':
            alphai = ((8 * U238 * (np.exp(lmbd238 * t1) - np.exp(lmbd238 * t2)))
                      + (7 * U235 * (np.exp(lmbd235 * t1) - np.exp(lmbd235 * t2)))
                      + (6 * Th232 * (np.exp(lmbd232 * t1) - np.exp(lmbd232 * t2))))
            
            anneal_d = alphai * self.damage
            anneal_d = np.sum(anneal_d, axis=0)
            
            El = 165
            D0l = 193188
            D0N17 = 0.0034
            EaN17 = 71
            Ba = 5.48E-19
            interconnection = 3
            SV = 1.669
            lint_lattice = 45920
        
            fa = 1 - np.exp(-Ba * anneal_d)
            DI = 1 - np.exp(-Ba * anneal_d * interconnection)
                        
            tort = (lint_lattice / (4.2 / (fa * SV) - 2.5)) ** 2
            Dtort = (1 / tort) * D0l * np.exp(-El / (R * T_mean))
            Dtorta2 = Dtort / (self.r * 1e-4 * (1 - DI))**2
            
            DN17 = D0N17 * np.exp(-EaN17 / (R * T_mean))
            DN17a2 = (DN17) / (self.r * 1e-4 * DI)**2
            
            diffusivities = (DI / DN17a2 + (1 - DI) / Dtorta2)**-1
            
            self.diffusivities = diffusivities * self.r**2
        
    def thermal_diffusion(self):
        """Diffusion (Farley et el., 2000) calculated as a function of
        temperature, with no rdaam. Runs faster, and may be appropriate for
        geologically rapid and recent cooling.
        """
        if self.mineral == 'apatite':
            Do = 50 #cm2/s
            Ea = 137.522 #kj/mol
        if self.mineral == 'zircon':
            Do = 0.46
            Ea = 169.0336
        
        R = 0.00831447
        T = self.T + 273.15
            
        d = (Do * np.exp(-Ea / (R * T))) * 1e8
        d = (d[:-1] + d[1:]) / 2 #average for dt
        
        self.diffusivities = d  

    def diffusion_solve(self):
        """Calculates a model cooling age given grain data, and a 
        time-temperature path.
        
        Returns
        -------
        Models age in million years
        """
        #these will be replaced by global variables
        nodes = self.nodes
        rstep = self.rstep
        x = self.x
        r = self.r

        ej238 = self.eject['ej238']
        ej235 = self.eject['ej235']
        ej232 = self.eject['ej232']

        t = self.t
        t1 = self.t1
        t2 = self.t2
        
        lmbd238 = self.lmbd238
        lmbd235 = self.lmbd235
        lmbd232 = self.lmbd232

        #production term A in atoms/g
        A = (8 * ej238[:,None] * (np.exp(lmbd238 * t1) - np.exp(lmbd238 * t2)) 
             + 7 * ej235[:,None] * (np.exp(lmbd235 * t1) - np.exp(lmbd235 * t2)) 
             + 6 * ej232[:,None] * (np.exp(lmbd232 * t1) - np.exp(lmbd232 * t2)))
        
        k = self.diffusivities

        def _1D_solver():
            """Solves the 1D diffusion eqn in sphere, following the approach of
            Ketcham, 2005, which uses the crank-nicolson method.
            """
            dt = (t1-t2) * 365.25 * 24 * 60 * 60
            dr = np.full((nodes,1),self.rstep)
            dr[0] = dr[0] / 2 # dr/2 at center
            
            beta = np.zeros((np.shape(A)))
            beta = beta + ((2 * dr**2) / (k * dt))
            
            ad0 = (-2 - beta)
            bd0 = (2 - beta)
            
            ad0[:,0] = - 3 - beta[:,0]
            bd0[:,0] = 3 - beta[:,0]

            kd1 = np.ones((len(beta)-1))
            
            u = np.full((nodes),100.0,dtype=np.float64)
            
            f = (A * beta * x[:,None])
            
            a = sparse.diags([kd1,ad0[:,0],kd1],[-1,0,1],format="csr")
            b = sparse.diags([-kd1,bd0[:,0],-kd1],[-1,0,1],format="csr")
            
            for i in range(len(t1)):      
                a.setdiag(ad0[:,i]) #update main diagonal for changing beta
                b.setdiag(bd0[:,i])
                u = spsolve(a, b.dot(u) - f[:,i])

            return u

        u = _1D_solver()

        integral = (u / x) * 4 * np.pi * x**2
        totalHe = integrate.romb(integral,dx=rstep)
        totalHe = totalHe / (4 * np.pi / 3)

        total238 = 8 * self.U238 * r**3
        total235 = 7 * self.U235 * r**3
        total232 = 6 * self.Th232 * r**3
        
        ft238 = self.eject['ft238']
        ft235 = self.eject['ft235']
        ft232 = self.eject['ft232']
        
        #iterative age calculation
        left_sum = (total238 * ft238
                    + total235 * ft235
                    + total232 * ft232
                    + totalHe
                    )
        
        lo_age = 0
        hi_age = t[0] * 1e6
        
        while (hi_age - lo_age) > 100:
            mid_age = (hi_age + lo_age) / 2
            
            mid_val = (total238 * ft238 * np.exp(lmbd238 * mid_age)
                       + total235 * ft235 * np.exp(lmbd235 * mid_age)
                       + total232 * ft232 * np.exp(lmbd232 * mid_age)
                       )
            
            if mid_val < left_sum:
                lo_age = mid_age
            else:
                hi_age = mid_age

        model_age = ((hi_age + lo_age) / 2) / 1e6

        return model_age
    
    def solve(self,t,T,k='rdaam'):
        """Main function to calculate He age for an inputted tT path.
        Diffusivity can be generated by RDAAM (Flowers et al., 2009), or simple
        diffusion (Farley et al., 2000) with no damage or annealing. The latter
        runs faster, but is inappropriate for most cooling histories.

        Inputs
        ------
        t:
            Time in million years
        T:
            Temperature in degrees C
        k:
            Diffusion model 'rdaam' (default) or 'simple'

        Returns
        -------
        Model He age [Ma]
        """
        startTime = runningtime.time()
        
        if np.max(T) < 80:
            print('\n Warning: max T should be >80 degrees C' +
                  ' to ensure model convergence')
        if (len(t)!=len(T)):
            print('\n Error: t and T arrays must be equal in size \n')
            sys.exit()
        if np.all(np.diff(t) > 0):
            print('\n Error: t values must decrease monotonically')
            sys.exit()

        self.T = T
        self.t = t
        self.t1 = self.t[:-1] * 1e6
        self.t2 = self.t[1:] * 1e6
                
        self.T_mean = (self.T[1:] + 273.15 + self.T[:-1] + 273.15) / 2

        if k == 'rdaam':
            self.rdaam_diffusivity()
            #self.t, self.T = self._get_dt()
            #self.rdaam_diffusivity()

        elif k == 'simple':
            self.thermal_diffusion()
        else:
            print('\n Error: k_model rdaam or diff must be selected 0 \n')
            sys.exit()
        
        He_age = self.diffusion_solve()
        print('\n .solve() ran in {:.2f}'.format(runningtime.time() - startTime))
        print('\n Helium age: {:.1f} Ma'.format(He_age))

        return He_age
