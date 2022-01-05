from __future__ import print_function, division
import numpy as np
import pandas as pd
from numpy.polynomial import hermite
from math import sqrt, pi

def smooth(x, sigma):
    '''smooth with a gaussian kernel'''


    kernel = np.exp(-.5 * np.arange(int(-4 * sigma), int(4 * sigma + 1), 1) ** 2 / sigma ** 2)
    kernel = kernel / sum(kernel)
    return np.convolve(x, kernel, 'same')

### Bayesian grid deconvolution
class DeltaDeconvolute:
    '''
    Deconvolute an underlying "signal" density from noisy measurements of data
    '''
    def __init__(self, noise_function, observed, min_support=1e-5, spacing=1, verbose=False):
        '''
        Parameters:
            noise_function: 1D array giving the value of the noise function at a set of points centered on
                the FIRST point.
                If points are not evenly spaced, the inference will not be correct. The inference
                assumes periodic boundary conditions, so use enough points to pad the observed data
            true_function: 1D array giving the density of the observed data at the same points.
            min_support (float): DeltaDeconvolute will only infer nonzero signal density where the observed
                data convoluted with a delta function is greater than this threshold
            spacing (int): evaluate the signal density every X many points. Must divide the points evenly
            verbose: get lots of information for debugging
        '''
        assert not any(noise_function < 0), 'Negative values in noise function'
        assert not any(observed < 0), 'Negative values in observed function'
        self.noise_function = noise_function/sum(noise_function)
        self.observed = observed/sum(observed)
        self.densities = []
        self.n = len(noise_function)
        assert not self.n%spacing, "spacing mismatched"
        assert self.n == len(observed), 'mismatched data!'
        self.evaluated_points = [True]*self.n
        for i in range(self.n):
            if i%spacing == 0:
                if sum(self.noise_function * self.observed) > min_support:
                    self.densities.append(
                     np.copy(self.noise_function)
                    )
                else:
                    self.evaluated_points[i] = False
            else:
                self.evaluated_points[i] = False
            self.noise_function = np.roll(self.noise_function, 1)
        self.lambdas = np.ones(len(self.densities), np.float)/len(self.densities)
        self.densities = np.stack(self.densities)
        self.q = np.ones(self.densities.shape, np.float)/len(self.densities)
        if verbose:
            print('n:', self.n)
            print('densities\n', self.densities)
            print('evaluated\n', self.evaluated_points)
            print()
    def fit(self, tol=1e-6, maxit=1000, verbose=False):
        self.cost = self.get_cost()
        for i in range(maxit):
            #E step: infer how much of each observation belongs to each density
            self.q[:] = np.outer(self.lambdas, self.observed) * self.densities / (
                np.sum(self.lambdas.reshape((-1, 1))*self.densities + 1e-16, axis=0)
                )
            if verbose:
                print('q:\n', self.q)
            #M step: infer size of each component
            self.lambdas = np.mean(self.q, axis=1)
            self.lambdas /= sum(self.lambdas)
            if verbose:
                print('lambda\n', self.lambdas)
            new_cost = self.get_cost()
            if verbose:
                print('cost:',  new_cost)
                print()
                last_change = (new_cost - self.cost)/self.cost
            self.cost = new_cost
            if 0 < last_change < tol:
                print("Converged at cost %f with %i iterations" %(new_cost, i))
                break
        if i >= maxit - 1:
            print("Reached max iterations %i with cost %f, last change %f" %(
                maxit, self.cost, last_change))

    def get_cost(self):
        return np.sum(np.square(np.convolve(self.noise_function, self.g, 'same') - self.observed))

    @ property
    def g(self):
        '''Get the inferred deconvoluted distribution. In general it will be shifted, often by half the
            the array length'''
        out = np.zeros(self.n)
        out[self.evaluated_points] = self.lambdas
        return out


### Differential deconvolution
def northdecomp(x, basis, maxit=1000, tol=1e-6):
    out = np.zeros((maxit, len(basis)))
    norm = sum(x**2)
    x = np.copy(x)/sqrt(sum(x**2))
    for i in range(maxit):
        for j, b in enumerate(basis):
            curr = sum(x*b)
            x -= curr*b
            out[i, j] = curr
        norm2 = sum(x**2)
        if (norm - norm2)/norm < tol:
            maxit = i
            break
        norm = norm2
    print('%i iterations, final norm: %f' %(maxit, norm))
    return np.sum(out, axis=0)

def differential_basis(x, n, step, dx=1, norms=None, smoothing=0, filt='linear'):
    x = np.copy(x)
    x /= sqrt(sum(x**2))
    new_norms = [1]
    out = [np.copy(x)]
    for i in range(n):
        x = differential(x, step, dx)
        if norms is None:
            new_norms.append(sqrt(sum(x**2)))
        else:
            new_norms.append(norms[i+1])
        x /= new_norms[-1]
        if filt == 'linear':
            x *= (n-i-1)*1.0/(n-1)
        if smoothing > 0:
            x = ut.smooth(x, smoothing)
        out.append(np.copy(x))
    return np.stack(out), new_norms

def hermite_basis(x, points, n):
    x = np.copy(x)
    out = []
    for i in range(n+1):
        out.append(
            hermite.hermval(points, [0]*(i)+[1]) *
            1.0/(sqrt(2**i)*sqrt_fact(i)) *
            x
        )
    return np.stack(out)


### Fourier deconvolution
def fourier(x, points, n):
    out = []
    norm = 1.0j*2*pi/(max(points) - min(points))
    for i in range(-n, n+1):
        out.append(np.sum(x*np.exp(i*norm*points)))
    return np.array(out)

def inverse_fourier(q, points):
    out = []
    norm = -1.0j*2*pi/(max(points) - min(points))
    coeffs = norm*np.arange(-len(q) // 2 + 1, len(q) // 2+1)
    for p in points:
        out.append(np.sum(q*np.exp(p*coeffs)))
    return out

def fourier_deconvolute(observed, noise, points, cutoff):
    f1 = fourier(observed, points, cutoff)
    f2 = fourier(noise, points, cutoff)
    factors = f1/f2
    triangle = 2 - (
        cutoff + 
        abs(np.arange(len(factors)) - .5*len(factors))
    ) / cutoff 
    triangle[triangle < 0] = 0
    return np.real(inverse_fourier(triangle*factors, points))
    


