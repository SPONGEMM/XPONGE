"""
This **module** implements the WHAM method
"""
from glob import glob
import numpy as np

class WHAM:
    """
        do the WHAM analysis

        :param window_edges: the edge of the windows
        :param temperature: the temperature
        :param weight: the weight of the bias
        :param references: the references of the bias
        :param period: the period of the CV
        :param step_limit: the maximum step to calculate the free energy
        :param diff_limit: stop iterations when the difference reaches the diff_limit 
    """
    def __init__(self, window_edges, temperature, weight, references, period=None, step_limit=3000, diff_limit=1e-4):
        self.diff_limit = diff_limit
        self.step_limit = step_limit
        self.beta = 4184 / 8.314 / temperature
        self.weight = weight
        self.references = references
        self.window_edges = window_edges
        self.cvs = None
        self.sampling_steps = None
        self.period = period

    def get_data_from_mdout(self, mdouts, cv_name):
        """
            get the CV information from the mdout files

            :param mdouts: the file name of mdouts
            :param cv_name: the name of the CV
        """
        from . import MdoutReader
        f = glob(mdouts)
        cvs = []
        sampling_steps = None
        for fi in f:
            mout = MdoutReader(fi)
            cvs.append(getattr(mout, cv_name))
            if sampling_steps is not None and sampling_steps != len(cvs[-1]):
                raise NotImplementedError("the reweighting for the simulations with different steps is not implemented")
            elif sampling_steps is None:
                sampling_steps = len(cvs[-1])
        self.sampling_steps = sampling_steps
        self.cvs = np.array(cvs)

    def bias(self, weight, x, ref):
        """the function to get the bias"""
        dx = x - ref
        if self.period is not None:
            dx -= np.floor(dx / self.period + 0.5) * self.period
        return weight*dx**2 

    def main(self):
        windows = (self.window_edges[1:] + self.window_edges[:-1]) / 2
        f = np.zeros(np.size(self.references))
        f_record = f.reshape(1, -1)
        for _ in range(self.step_limit):
            bias = self.bias(self.weight, self.cvs.reshape(1, -1, self.sampling_steps), self.references.reshape(-1, 1, 1))
            numerator = np.exp(-self.beta*bias)
            denominator = np.sum(self.sampling_steps * np.exp(self.beta * f).reshape(-1,1,1) * numerator, axis = 0)
            f = -np.log(np.sum(numerator/denominator, axis = (1,2)))/self.beta
            f_record = np.vstack([f_record, f.reshape(1,-1)])
            if np.abs(np.max(f_record[-1,:] - f_record[-2,:])) < self.diff_limit:
                break

        prob = np.zeros_like(windows)
        for i in range(len(prob)):
            count = len(np.where(( self.cvs >= self.window_edges[i] ) & ( self.cvs < self.window_edges[i+1] ))[0])
            if count:
                bias = self.bias(self.weight, windows[i], self.references)
                prob[i] = count / np.sum(self.sampling_steps * np.exp(-self.beta * (bias - f)))
            else:
                prob[i] = 0

        prob = prob/np.sum(prob)
        free_energy = -np.log(prob)/self.beta
        free_energy -= np.min(free_energy)
        return windows, free_energy, f_record

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import random
    kT = 1.2
    MCsteps = 5000
    weight = 200
    sigma = 0.11

    G = 10
    U = lambda x: G*(x-1)**2*(x+1)**2
    lower_boundary = -1.5
    upper_boundary = 1.5

    num_windows = 20
    test_boundaries = np.linspace(lower_boundary, upper_boundary, num = num_windows+1)
    X_reference = test_boundaries[:-1] + 0.5*np.mean(np.diff(test_boundaries))
    X_record = np.zeros([num_windows, MCsteps])
    E_record = np.zeros([num_windows, MCsteps])
    X_current = X_reference[0]
    def MCMC(kT,
            MCsteps,
            X_reference,
            lower_boundary,
            upper_boundary,
            weight,
            U,
            X_current,
            sigma):
        beta = 1/kT
        U_biased = lambda x: U(x) + weight*(x - X_reference)**2
        E_current = U_biased(X_current)
        X_record = np.zeros(MCsteps)
        E_record = np.zeros(MCsteps)   
        X_record[0] = X_current
        E_record[0] = E_current    
        for i in range(1,MCsteps,1):   
            X_new = X_current + sigma*np.random.normal(0,1,1)  
            if X_new < lower_boundary or X_new > upper_boundary:     
                X_record[i] = X_current
                E_record[i] = E_current
                continue
            else:   
                E_new = U_biased(X_new)
                dE = E_new - E_current
                P = np.exp(-beta*dE)
                u = random.uniform(0, 1)        
                if u <= min(1,P):           
                    E_current = E_new
                    X_current = X_new          
                X_record[i] = X_current
                E_record[i] = E_current       
        return [X_record, E_record]
    count = 0
    for i in range(0,num_windows,1): 
        count = count + 1
        print(count)
        [ X_record[i,:], E_record[i,:] ] = MCMC(kT,
                                                MCsteps,
                                                X_reference[i],
                                                lower_boundary,
                                                upper_boundary,
                                                weight,
                                                U,
                                                X_current,
                                               sigma)  
        if i < num_windows-1:
            sampled_states_so_far = (np.asanyarray(X_record)).flatten()
            idx = (np.abs(sampled_states_so_far - X_reference[i+1])).argmin()
            X_current = sampled_states_so_far[idx]      
        else:        
            continue

    w = WHAM(np.linspace(lower_boundary, upper_boundary, 101), 600, 200, X_reference)
    w.cvs = X_record
    w.sampling_steps = MCsteps
    x, y, f_record = w.main()

    plt.plot(x, y)

    positions = np.linspace(-np.pi, np.pi, 150)
    G = 10
    U = lambda x: w.beta*G*(x-1)**2*(x+1)**2
    x = np.linspace(lower_boundary,upper_boundary,10**3)
    plt.plot(x,U(x))

    plt.show()
    plt.plot(f_record)
    plt.show()
