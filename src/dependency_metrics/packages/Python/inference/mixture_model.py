from __future__ import division, print_function
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import gaussian_kde, norm
from scipy.interpolate import interp1d


def smooth(x, sigma):
    '''smooth with a gaussian kernel'''


    kernel = np.exp(-.5 * np.arange(int(-4 * sigma), int(4 * sigma + 1), 1) ** 2 / sigma ** 2)
    kernel = kernel / sum(kernel)
    return np.convolve(x, kernel, 'same')

### Infer class probabilities
class MixFitOneUnknown:
    '''
    fit a 1D mixture model with a set of known functions (not necessarily gaussian)
    and an optional unknown gaussian component
    '''
    def __init__(self, densities, data, initial_lambdas=None, fit_gaussian=True,
            initial_mu=None, initial_sigma=None):
        '''
        Parameters:
            densities: iterable of normalized functions R^N -> P^N
            data: 1D array of points
            initial_lambdas: iterable with same length as densities giving initial guesses for size of each
                component. Must sum to between 0 and 1. Must sum to 1 if fit_gaussian=False
            fit_gaussian: Bool designating whether to include gaussian
            lambda_lock: if True, model will not update the mixing proportions of the components
        '''
        if initial_lambdas is None:
            initial_lambdas = [1.0/(len(densities) + int(fit_gaussian)) for d in densities]
        self.fit_gaussian = fit_gaussian
        if sum(pd.isnull(densities)) > 0:
            e = 'Error: %i null values in data\n%r' %(sum(pd.isnull(densities)), data)
            raise ValueError(e)
        if fit_gaussian:
            assert 0 < sum(initial_lambdas) < 1, "invalid mixing sizes (initial_lambdas)"
            self.lambdas = np.concatenate([ initial_lambdas, [1-sum(initial_lambdas)] ])
        else:
            self.lambdas = np.array(initial_lambdas)
        self.n = len(initial_lambdas)+1
        self.densities = list(densities)
        self.data = np.array(data)
        if fit_gaussian:
            self.densities.append(self.gauss)
            if initial_mu is None:
                self.mu = np.mean(data)
            else:
                self.mu = initial_mu
            if initial_sigma is None:
                self.sigma = np.std(data)
            else:
                self.sigma = initial_sigma
            self.sigma = np.std(data)
        else:
            assert sum(initial_lambdas) == 1, 'Invalid mixing sizes if not fitting gaussian'
        self.q = 1.0/(len(self.densities)) * np.ones((len(self.densities), len(self.data)))
        self.p = np.stack([d(self.data) for d in densities])
    
    def gauss(self, x):
        return norm.pdf(x, loc=self.mu, scale=self.sigma)
    
    def fit(self, tol=1e-7, maxit=1000, debug=False, lambda_lock=False):
        self.cost = self.get_cost()
        for i in range(maxit):

            #E step
            for k in range(self.n-1):
                self.q[k, :] = self.lambdas[k]*self.p[k, :]
            if self.fit_gaussian:
                self.q[self.n-1, :] = self.lambdas[self.n-1]*self.gauss(self.data)
            if any(np.sum(self.q, axis=0) == 0):
                loc = np.argwhere(np.sum(self.q, axis=0) == 0).ravel()
                bad_points = self.data[loc]
                e = 'All component densities invalid for indices %r\ndata there: %r\ndensities there: %r)' %(
                    loc, bad_points, [d(bad_points) for d in self.densities])
                raise ValueError(e)
            self.q /= np.sum(self.q, axis=0)
            if debug:
                assignment = np.argmax(self.q, axis=0)
                for i in range(len(self.lambdas)):
                    vals = self.data[assignment == i]
                    print('post-normalize: component %i min %f max %f' %(i, min(vals), max(vals)))

            #M step
            if self.fit_gaussian:
                self.mu = np.sum(self.q[-1]*self.data)/np.sum(self.q[-1])
                if debug: print('mu %f' %self.mu)
                self.sigma = np.sqrt(np.sum(
                    self.q[-1]*np.square(self.data-self.mu)
                    ) / np.sum(self.q[-1]))
                if debug: print('sigma %f' %self.sigma)
            if not lambda_lock:
                self.lambdas = np.mean(self.q, axis=1)
                self.lambdas /= sum(self.lambdas)
            if debug: print('lambdas %r' %self.lambdas)
            new_cost = self.get_cost()
            if debug: print('new cost %f' %new_cost )
            last_change = (new_cost - self.cost)/self.cost
            self.cost = new_cost
            if 0 < -last_change < tol:
                print("Converged at cost %f with %i iterations" %(new_cost, i))
                break
        if i >= maxit - 1:
            print("Reached max iterations %i with cost %f, last change %f" %(
                maxit, self.cost, last_change))
        
    def get_cost(self):
        out = -np.mean(
            np.log(np.sum([q*d(self.data)+1e-32
            for q, d in zip(self.q, self.densities)
        ])))
        if np.isnan(out):
            e = "Null cost. %i nulls and %i negatives in q values, %i nulls and %i negatives in densities\
            \n%r\n%r" %(
                np.sum(np.isnan(self.q)),
                np.sum(self.q < 0),
                np.sum(np.isnan(np.stack([d(self.data) for d in self.densities]))),
                np.sum(np.stack([d(self.data) for d in self.densities]) < 0),
                np.argwhere(np.isnan(self.q)),
                self.q.shape
                )
            raise ValueError(e)
        return out
    
    def full_density(self, points):
        out = 0*points
        for d,l in zip(self.densities, self.lambdas):
            out += l*d(points)
        return out
    
    def get_assignments(self, component):
        '''get probability that each data point was generated by the given component (int)'''
        return self.q[component]

    def component_probability(self, points, component):
        '''get probability that data at specified points belongs to the given component (int)'''
        return self.lambdas[component] * self.densities[component](points) / sum([
            l*d(points) for l, d in zip(self.lambdas, self.densities)
            ])


class TailKernel:
    '''Used to set points at left and right extremes to fixed values'''
    def __init__(self,
                lower_bound=None, lower_value=None, upper_bound=None, upper_value=None):
        self.lower_bound = lower_bound
        self.lower_value = lower_value
        self.upper_bound = upper_bound
        self.upper_value = upper_value

    def apply(self, x, y):
        if self.lower_bound is not None:
            y[x < self.lower_bound] = self.lower_value
        if self.upper_bound is not None:
            y[x > self.upper_bound] =self.upper_value
                

def probability_2class(component0_points, component1_points, all_points,
               smoothing='silverman', p_smoothing=.15,
               right_kernel_threshold=None, left_kernel_threshold=None, probability_kernel=None,
               maxit=500, lambda_lock=False, **kwargs):
    '''
    Estimates the distributions of component0_points and component1_points using a gaussian kernel,
    then assigns each of all_points a probability of belonging to the component distribution. Note that
    this is NOT a p-value: P(component1) = 1 - P(component0)
    '''
    #estimate density
    estimates = [
            gaussian_kde(component0_points, bw_method=smoothing),
            gaussian_kde(component1_points, bw_method=smoothing)
        ]
    points = np.arange(min(all_points) - 4*p_smoothing, max(all_points) + 4*p_smoothing, .01)
    estimates = [e(points) for e in estimates]
    
    #this step copes with the fact that scipy's gaussian KDE often decays to true 0 in the far tails, leading to
    #undefined behavior in the mixture model
    if right_kernel_threshold is not None:
        estimates[0][np.logical_and(points > right_kernel_threshold, estimates[0] <1e-16)] = 1e-16
    if left_kernel_threshold is not None:
        estimates[1][np.logical_and(points < left_kernel_threshold, estimates[1] <1e-16)] = 1e-16

    #create density functions using interpolation (these are faster than calling KDE on points)
    densities = [interp1d(points, e) for e in estimates]
    
    #infer probabilities
    fitter = MixFitOneUnknown(densities, all_points, **kwargs)
    fitter.fit(maxit=maxit, lambda_lock=lambda_lock)

    #generate smoothed probability function
    p = fitter.component_probability(points, 1)
    if probability_kernel is not None:
        probability_kernel.apply(points, p)
    p = smooth(p, int(100*p_smoothing))
    d = interp1d(points, p)

    #return probabilities
    out = d(all_points)
    out[out > 1] =  1
    out[out < 0] = 0
    return out

def assign_probabilities(melted, line, lambdas, smoothing=.05, p_smoothing=.15, **kwargs):
    '''
    Given a melted matrix of CERES scores, assign each score a probability of belonging to the
    dependent class.
    args:
        melted: melted matrix of CERES scores. Must have the following columns:
            'lines': cell line names
            'genes': gene names
            'score': CERES score
            'unexpressed': Boolean, true if the cell line does not express that gene
            'essential': Boolean, true if the gene is believed to be a dependency in all lines
            'nonessential': Boolean, true if the gene is believed to be a dependency in no lines
            'p_dependent': Probability of the gene being a dependency for that cell line. This will be
                    overwritten in the specified cell line, and can be empty.
        line: name of line to assign probabilities in, or 'all'
        smoothing: width of gaussian kernel for estimating the distributions of the components
        p_smoothing: width of smoothing kernel for the probability values
        additional arguments passed to probability_2class
    '''           
    if line == 'all':
        subset = melted
    else:
        subset = melted[melted.lines == line]
    if subset.unexpressed.sum() > 0:
        component0_points = subset.score[subset.unexpressed]
    else:
        component0_points = subset.score[subset.nonessential]
    component1_points = subset.score[subset.essential]
    probability_kernel = TailKernel(lower_bound=-1.5, lower_value=1, upper_bound=.25, upper_value=0)
    melted.loc[subset.index, 'p_dependent'] = probability_2class(
        component0_points, component1_points, subset.score, smoothing=smoothing, p_smoothing=p_smoothing, 
        right_kernel_threshold=0, left_kernel_threshold=-.8, probability_kernel=probability_kernel,
        initial_lambdas=lambdas, **kwargs)      
        

def test():
    '''Assume we have some mixture where we know two components are T-distributions with 4 and 5 degrees
    of freedom, centered at -2 and 0 respectively. There's also a third component which we'll treat as
    a gaussian and allow the mixture model to fit.
    '''
    #fake data:
    print('testing with fake T-distributed data')
    kwargs1 = {'df': 4, 'loc': -2, 'scale': .5}
    kwargs2 = {'df': 5, 'loc': 0, 'scale': .25}
    kwargs3 = {'df': 15, 'loc': -1, 'scale': 1}
    test_data = np.concatenate([
        stats.t.rvs(size=10000, **kwargs1),
        stats.t.rvs(size=60000, **kwargs2),
        stats.t.rvs(size=30000, **kwargs3)
        ])
    np.random.shuffle(test_data)
    
    def component_0(x):
        return stats.t.pdf(x, **kwargs1)
    def component_1(x):
        return stats.t.pdf(x, **kwargs2)
    
    #Initial guesses for how much each of the data belongs to each of the known density functions:
    lambdas = [.2, .3]
    
    #fit:
    mixture = MixFitOneUnknown(densities=[component_0, component_1], data=test_data,
                               initial_lambdas=lambdas)
    print('training mixture model')
    mixture.fit(maxit=1000, tol=1e-8)
    print('mean: %f (expected %f), std: %f (expected %f), lambdas: %r (expected [.1, .6, .3])' %(
    mixture.mu, kwargs3['loc'], mixture.sigma, kwargs3['scale'], mixture.lambdas))


