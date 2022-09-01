#%%
import numpy as np
from inspect import signature
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
from scipy.stats import chi2
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#%%
''' Utillity functions '''
from scipy.special import gamma
from scipy.special import digamma
from scipy.special import polygamma


# Students-t distribution

def t_dist(x, a, mu, sigma, nu, b):
    num = gamma((nu + 1) / 2)
    denom = np.sqrt(nu * np.pi) * sigma * gamma(nu / 2)
    exp = 1 + ((x - mu) / sigma) ** 2 / nu
    pow = - (nu + 1) / 2
    return a * num / denom * exp ** pow + b

# parameter gradient
def t_dist_jac(x, a, mu, sigma, nu, b):
    exp = 1 + ((x - mu) / sigma) ** 2 / nu
    exp_deriv = (nu + 1) * (x - mu) / sigma ** 2 / nu
    nu_dev_1 = (nu + 1) * (x - mu) ** 2 / (2 * sigma ** 2 * nu ** 2)
    nu_dev_2 = (digamma((nu + 1) / 2) - digamma(nu / 2)) / 2

    a_deriv = t_dist(x, a, mu, sigma, nu, 0) / a

    mu_deriv = exp_deriv * exp ** (-1) * t_dist(x, a, mu, sigma, nu, 0)

    sigma_deriv = ((x-mu)**2 - sigma**2) / sigma ** 3 * exp ** (-1) * t_dist(x, a, mu, sigma, nu, 0)

    nu_deriv = (nu_dev_1 - np.log(exp) + nu_dev_2) * t_dist(x, a, mu, sigma, nu, 0)

    if isinstance(x, np.ndarray) != 1:
        b_deriv = 1
    else:
        b_deriv = np.ones(len(x))
    return np.array([a_deriv, mu_deriv, sigma_deriv, nu_deriv, b_deriv])

# Hessian
def t_dist_hess_matrix(x, a, mu, sigma, nu, b):
    f0 = t_dist(x, a, mu, sigma, nu, 0)
    exp = 1 + ((x - mu) / sigma) ** 2 / nu
    exp_deriv = (nu + 1) * (x - mu) / sigma ** 2 / nu
    denom = ((x - mu) ** 2 + sigma ** 2 * nu)
    jacob = t_dist_jac(x, a, mu, sigma, nu, b)

    grad_deriv_a = np.insert(np.insert(jacob[1:-1] / a, 0, 0), 4, 0)

    grad_deriv_mu = np.insert(np.insert(jacob[1:-1] * exp_deriv / exp +
                              f0 * np.array([(nu+1) * ((x-mu)**2-sigma**2*nu) / denom**2,
                                             exp_deriv / exp * 2 * sigma * nu / denom,
                                             ((x-mu)**2 - sigma**2) * (x - mu) / denom**2]),
                                             0, jacob[1] / a), 4, 0)

    grad_deriv_sig = np.insert(np.insert(jacob[1:-1] * ((x-mu)**2 - sigma**2) / sigma ** 3 / exp +
                                f0 * np.array([2 * (nu+1)*(mu-x) * sigma * nu / denom**2,
                                               ((nu*sigma)**2 - nu*(x-mu)**4/sigma**2 - nu*(3*nu+1)*(x-mu)**2) / denom**2,
                                               - (x - mu) ** 2 *((x - mu) ** 2 + sigma**2) / sigma / denom**2]),
                                         0, jacob[2] / a), 4, 0)

    grad_deriv_nu = np.insert(np.insert(jacob[1:-1] * jacob[3] / f0 +
                              f0 * np.array([2 * (x - mu) / denom - (nu + 1) * (x - mu) / (sigma * nu) ** 2,
                                             2 * (x - mu) ** 2 / denom / sigma - (nu + 1) * (x - mu) ** 2 / sigma ** 3 / nu ** 2,
                                            (1/nu - (x-mu)**2/2/(nu*sigma)**2 -
                                             sigma**2/denom -
                                             (x-mu)**2/nu**3/sigma**3) +
                                             (polygamma(1, nu / 2) - polygamma(1,(nu + 1) / 2)) / 4]),
                                             0, jacob[3] / a), 4, 0)

    grad_deriv_b = np.array([0, 0, 0, 0, 1])

    return np.array([grad_deriv_a, grad_deriv_mu,
                     grad_deriv_sig, grad_deriv_nu,
                     grad_deriv_b])

def stack(A, B):
    Ae = 1
    Be = 1
    if isinstance(A[0, 0], np.ndarray):
        Ae = len(A[0, 0])
    if isinstance(B[0, 0], np.ndarray):
        Be = len(B[0, 0])
    C = np.empty((A.shape[0], A.shape[1], Ae + Be))
    if Ae == 1 and Be == 1:
        for i in range(len(A)):
            for j in range(len(A[i])):
                C[i, j] = np.array([A[i, j], B[i, j]])
    elif Ae != 1 and Be==1:
        for i in range(len(A)):
            for j in range(len(A[i])):
                C[i, j] = np.append(A[i, j], B[i, j])
    else:
        for i in range(len(A)):
            for j in range(len(A[i])):
                C[i, j] = np.concatenate((A[i, j], B[i, j]))
    return C

def t_dist_hess(x, a, mu, sigma, nu, b):
    if isinstance(x, np.ndarray):
        ret = t_dist_hess_matrix(x[0], a, mu, sigma, nu, b)
        for i in range(1, len(x)):
            ret = stack(ret, t_dist_hess_matrix(x[i], a, mu, sigma, nu, b))
        return np.array(ret)
    else:
        return t_dist_hess_matrix(x, a, mu, sigma, nu, b)

# Gaussian distribution

def gauss_dist(x, a, mu, sigma, b):
    norm = 1 / np.sqrt(2*np.pi) / sigma
    pow = - (x - mu) ** 2 / 2 / sigma**2
    return a * norm * np.exp(pow) + b

# parameter gradient
def gauss_dist_jac(x, a, mu, sigma, b):
    f0 = gauss_dist(x, a, mu, sigma, 0)
    a_deriv = f0 / a
    mu_deriv = f0 * (x - mu) / sigma**2
    sigma_deriv = f0 * (x - mu)**2 / sigma**3
    if isinstance(x, np.ndarray) != 1:
        b_deriv = 1
    else:
        b_deriv = np.ones(len(x))
    return np.array([a_deriv, mu_deriv, sigma_deriv, b_deriv])

# Hessian
def gauss_dist_hess_matrix(x, a, mu, sigma, b):
    f0 = gauss_dist(x, a, mu, sigma, 0)
    jacob = gauss_dist_jac(x, a, mu, sigma, b)

    grad_deriv_a = np.insert(np.insert(jacob[1:-1] / a, 0, 0), 3, 0)

    grad_deriv_mu = np.insert(np.insert(jacob[1:-1] * (x-mu) / sigma**2 -
                                        f0 * np.array([1/sigma**2, 2*(x-mu)/sigma**3])
                              , 0, jacob[1] / a), 3, 0)

    grad_deriv_sig = np.insert(np.insert(jacob[1:-1] * (x-mu)**2 / sigma**3 -
                                         f0 * np.array([2*(x-mu)/sigma**3, 3*(x-mu)**2/sigma**4])
                               , 0, jacob[1] / a), 3, 0)

    grad_deriv_b = np.array([0, 0, 0, 1])
    return np.array([grad_deriv_a, grad_deriv_mu,
                     grad_deriv_sig, grad_deriv_b])

def gauss_dist_hess(x, a, mu, sigma, b):
    if isinstance(x, np.ndarray):
        ret = gauss_dist_hess_matrix(x[0], a, mu, sigma, b)
        for i in range(1, len(x)):
            ret = stack(ret, gauss_dist_hess_matrix(x[i], a, mu, sigma, b))
        return np.array(ret)
    else:
        return gauss_dist_hess_matrix(x, a, mu, sigma, b)


#%%

def masked_matrix(arr,mask):
    deleted = 0
    for i in range(len(mask)):
        if not mask[i]:
            arr = np.delete(arr, i-deleted, 0)  # axis: 0 is row, 1 is col
            arr = np.delete(arr, i-deleted, 1)  # axis: 0 is row, 1 is col
            deleted += 1
    return arr

class FittingRoutine:

    """
    This class is a wrapper for scipy functions used to do curvefitting.

    Parameters
    ----------
    function : callable
        The objective function to be minimized.
            ``f(x, *params) -> float``
        where ``x`` is a 1-D array with shape (n,) and ``params``
        is a tuple of the fixed parameters needed to completely
        specify the function.

    x : ndarray, shape (n,)
        The independent variable where the data is measured of size (n,).

    y : ndarray, shape (n,)
        The dependent data, a length n array - nominally ``f(xdata, ...)``.

    error_y: ndarray, shape (n,), shape (2,n), optional
        Determines the uncertainty in `y`. Can either be
        of shape (n,) for symmetric uncertainty or (2,n) for
        assymetric uncertainty.

    error_x : ndarray, shape (n,), optional
        Determines the uncertainty in `x` of size (n,).

    P0 :  ndarray, shape (m,), shape (m,2), optional
        Initial guess. Array of real elements of size (m,),
        where ``m`` is the number of fixed parameters given to function.
        If ``method`` given as 'diff_evol' the search area should be given
        as array of real elements of size (m,2)

    covy : ndarray, shape (n,n), shape (2,n,n), optional
        Covariance matrix between ``n`` ``y`` points, either
        symmetric shape (n,n) or assymetric shape (2,n,n)

    cov_x : ndarray, shape (n,n), optional
        Covariance matrix between ``n`` ``x`` points of shape (n,n)

    method : str, optional
        Type of solver.  Should be one of
            - 'BFGS'        :ref:`<optimize.minimize-bfgs>`
            - 'Newton-CG'   :ref:`<optimize.minimize-newtoncg>`
            - 'dogleg'      :ref:`<optimize.minimize-dogleg>`
            - 'trust-ncg'   :ref:`<optimize.minimize-trustncg>`
            - 'trust-exact' :ref:`<optimize.minimize-trustexact>`
            - 'trust-krylov' :ref:`<optimize.minimize-trustkrylov>`
            - 'diff_evol' :ref:`<optimize.differential_evolution>`

        Standard is 'BFGS'. For any other method 'jac' should be specified

    mask : ndarray, shape (n,), optional
        Array of bolean values choosing points 'x,y' to use for 'True'
        and which not to use for 'False' of shape (n,)

    jac : callable
        If it is a callable, it should be a function that returns the gradient
        vector:
            ``jac(x, *args) -> array_like, shape (m,)``
        of the ``m`` fixed parameters for ``function``,
        where ``x`` is an array with shape (n,) and ``args`` is a tuple with
        the fixed parameters.

    Returns
    -------
    res : Result
        FitResult : Scipy OptimizeResult
            The optimization result represented as a `OptimizeResult` object.
            Important attributes are: ``x`` the solution array, ``success`` a
            Boolean flag indicating if the optimizer exited successfully and
            ``message`` which describes the cause of the termination.
        params : Fixed funtion parameters for minimized value of Chi2
        df : Degrees of freedom of minimization problem
        Chi2 :  Minimized value of Chi2
        Cov :  Covariance matrix of parameters for minimized value of Chi2
        Pval :  P-value (survival function of Chi2(Chi2,df)) for minimized value of Chi2
        NormRes : Normalized residuals for minimized value of Chi2
    Notes
    -----
    Fitting is achieved by minimizing:
    .. math::
        \chi^2 = \left(\mathbf{y} - f(\mathbf{x},*params))\left[\Sigma_y^c+\Sigma_x^cJ(\mathbf{x})\right]^{-1}(\mathbf{y} - f(\mathbf{x},*params))\right|

    Examples
    --------
    >>> import numpy as np
    >>> from FittingRoutine import FittingRoutine
    >>> def func(x, a, b):
    ...     return a/np.sqrt(2*np.pi) * np.exp(- x**2/2) + b
    >>> def func_deriv(x, a, b):
    ...     a_deriv = 1/np.sqrt(2*np.pi) * np.exp(- x**2/2)
    ...     b_deriv = np.ones(len(x))
    ...     return np.array([a_deriv,b_deriv])
    Define the data to be fit with some noise from covariance matrix:
    >>> N = 100
    >>> xdata = np.linspace(-3, 3, N)
    >>> rng = np.random.default_rng()
    >>> y_noise = rng.normal(0,1,size=xdata.size)
    >>> ycov = np.ones((N,N))*0.8
    >>> np.fill_diagonal(ycov, 1)
    >>> mean = func(xdata, 20, 0)+y_noise
    >>> ydata = np.random.multivariate_normal(mean=mean, cov=ycov, size=1)[0]
    Fit for the parameters a, b of the function `func`:
    >>> Fit = FittingRoutine(func, xdata, ydata, covy=ycov, jac=func_deriv)
    >>> print(Fit.params)
    [19.66767691 -0.09370664]
    >>> print(Fit.Cov)
    [[ 0.10524428 -0.01737613]
    [-0.01737613  0.00388857]]
    >>> print(Fit.Pval)
    0.6628323592625245
    """

    def __init__(self, function, x, y, error_y = np.array([None]), error_x = np.array([None]), P0 = np.array([None]),
                 covy=np.array([None]), covx=np.array([None]), method = None, mask = np.array([None]), jac=None, hess=None, **kwargs):

        '''Check if function is a function'''
        if callable(function) != 1:
            raise SyntaxError('First argument must be a callable function')
        else:
            pass
        if jac != None:
            if callable(jac) != 1:
                raise SyntaxError('Jac argument must be a callable function')
            else:
                pass

        '''Check if function is compatible'''
        if len(signature(function).parameters) <= 1:
            raise SyntaxError('Function has to few arguments')
        else:
            pass
        if jac != None:
            if len(signature(jac).parameters) != len(signature(function).parameters):
                raise SyntaxError('Function and jacobian has different number of parameters')
            else:
                pass
        if hess != None:
            if len(signature(hess).parameters) != len(signature(function).parameters):
                raise SyntaxError('Function and hessian has different number of parameters')
            elif jac == None:
                print('Hessian will not be used when jacobian is not given')
            else:
                pass

        '''Check if data has the right format'''
        if isinstance(x, (np.ndarray)) != 1:
            raise SyntaxError('x has to be a numpy array')
        else:
            pass
        if isinstance(y, (np.ndarray)) != 1:
            raise SyntaxError('y has to be a numpy array')
        else:
            pass
        if isinstance(error_x, (np.ndarray)) != 1:
            raise SyntaxError('error_x has to be a numpy array')
        else:
            self.error_x = error_x
        
        if isinstance(error_y, (list, tuple, np.ndarray)) != 1:
            raise SyntaxError('error_y has to be a numpy array')
        else:
            self.error_y = error_y

        if isinstance(covy, (list, tuple, np.ndarray)) != 1:
            raise SyntaxError('covy has to be a numpy array')
        else:
            self.cov_y = covy
            self.cov_y_plot = covy

        if isinstance(covx, (list, tuple, np.ndarray)) != 1:
            raise SyntaxError('covx has to be a numpy array')
        else:
            self.cov_x = covx
            self.cov_x_plot = covx

        if isinstance(mask, np.ndarray) != 1:
            raise SyntaxError('mask has to be a numpy array')
        else:
            self.mask = np.array([True for i in range(len(x))])

        '''  Check compatibility of x and y'''
        if (len(x) - len(y)) != 0:
            raise SyntaxError('Length of data does not match. Data of length ' + str(
                len(x)) + ' is not compatible with data of length ' + str(len(y)))

        ''' Check compatibility of mask'''
        if mask.any() != None:
            if len(mask) - len(x) != 0:
                raise SyntaxError(f'Length of mask {len(mask)} does not match that of data {len(x)}')
            else:
                self.mask = mask

        '''Check if data is compatible'''
        if covx.any() != None or covy.any() != None:
            # Use covariance matrix instead of normal error.
            if (covx.any() != None and covy.any() == None):
                raise NotImplementedError("Only errors on x not implemented yet")
                if len(covx.shape) != 2:
                    raise SyntaxError(f'Covariance matrix covx has to be (n,n) matrix not {covx.shape}')
                elif (covx.shape[0] - len(x)) != 0 or \
                        (covx.shape[1] - len(x)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(
                        len(x)) + ' is not compatible with data of length ' + str(covx.shape))
                else:
                    self.cov_x = masked_matrix(self.cov_x,self.mask)

            elif (covx.any() != None and covy.any() != None):
                if len(covx.shape) != 2:
                    raise SyntaxError(f'Covariance matrix covx has to be (n,n) matrix not {covx.shape}')
                elif (covx.shape[0] - len(x)) != 0 or \
                        (covx.shape[1] - len(x)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(
                        len(x)) + ' is not compatible with data of length ' + str(covx.shape))
                else:
                    self.cov_x = masked_matrix(self.cov_x,self.mask)
                if len(covy.shape) > 3 or len(covy.shape) < 2:
                    raise SyntaxError(f'Covariance matrix covy has to be (n,n) or (n,n,2) matrix not {covx.shape}')
                elif len(covy.shape) == 3:
                    if (covy.shape[1] - len(y)) != 0 or \
                            (covy.shape[2] - len(y)) != 0:
                        raise SyntaxError('Length of data does not match. Data of length ' + str(
                            len(y)) + ' is not compatible with data of length ' + str(covy.shape))
                    elif covy.shape[0] != 2:
                        raise SyntaxError(f'Covariance matrix covy has to be (n,n) or (2,n,n) matrix not {covx.shape}')
                    else:
                        self.cov_y = np.array([masked_matrix(self.cov_y[0],self.mask),masked_matrix(self.cov_y[1],self.mask)])
                        self.inv_cov_y = np.array(
                            [np.abs(np.linalg.inv(self.cov_y[0])), np.abs(np.linalg.inv(self.cov_y[1]))])
                elif (covy.shape[0] - len(y)) != 0 or \
                        (covy.shape[0] - len(y)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(
                        len(y)) + ' is not compatible with data of length ' + str(covy.shape))
                else:
                    self.cov_y = masked_matrix(self.cov_y,self.mask)
                    self.inv_cov_y = np.abs(np.linalg.inv(self.cov_y))

            elif (covx.any() == None and covy.any() != None):
                if len(covy.shape) > 3 or len(covy.shape) < 2:
                    raise SyntaxError(f'Covariance matrix covy has to be (n,n) or (n,n,2) matrix not {covx.shape}')
                elif len(covy.shape) == 3:
                    if (covy.shape[1] - len(y)) != 0 or \
                         (covy.shape[2] - len(y)) != 0:
                        raise SyntaxError('Length of data does not match. Data of length ' + str(
                            len(y)) + ' is not compatible with data of length ' + str(covy.shape))
                    elif covy.shape[0] != 2:
                        raise SyntaxError(f'Covariance matrix covy has to be (n,n) or (2,n,n) matrix not {covx.shape}')
                    else:
                        self.cov_y = np.array([masked_matrix(self.cov_y[0],self.mask),masked_matrix(self.cov_y[1],self.mask)])
                        self.inv_cov_y = np.array([np.abs(np.linalg.inv(self.cov_y[0])), np.abs(np.linalg.inv(self.cov_y[1]))])
                elif (covy.shape[0] - len(y)) != 0 or \
                        (covy.shape[0] - len(y)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(
                        len(y)) + ' is not compatible with data of length ' + str(covy.shape))
                else:
                    self.cov_y = masked_matrix(self.cov_y,self.mask)
                    self.inv_cov_y = np.abs(np.linalg.inv(self.cov_y))

        else:
            if (error_x.any() != None and error_y.any() == None):
                raise NotImplementedError("Only errors on x not implemented yet")
                if (len(x)-len(error_x)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(len(x)) + ' is not compatible with data of length ' + str(len(error_x)))
                else:
                    self.cov_x_plot = np.diag(self.error_x ** 2)
                    self.error_x = self.error_x[self.mask]
                    self.cov_x = np.diag(self.error_x**2)


            elif (error_x.any() == None and error_y.any() != None):
                if len(error_y.shape) > 1:
                    if error_y.shape[0] == 2 and isinstance(error_y.shape[1], int):
                        if (len(x)-len(error_y[0])) != 0 or (len(x)-len(error_y[1])) != 0:
                            raise SyntaxError('Length of data does not match. Data of length ' + str(len(x)) + ' is not compatible with data of length ' + str(error_y.shape))
                        else:
                            self.cov_y_plot = np.array([np.diag(self.error_y[0] ** 2), np.diag(self.error_y[1] ** 2)])
                            self.error_y = self.error_y[:,self.mask]
                            self.cov_y = np.array([np.diag(self.error_y[0]**2),np.diag(self.error_y[1]**2)])
                            self.inv_cov_y = np.array([np.abs(np.linalg.inv(self.cov_y[0])), np.abs(np.linalg.inv(self.cov_y[1]))])
                    else:
                        raise SyntaxError('Wrong shape for error_y ' + str(error_y.shape) + ' must be of shape (2,n) or (n,)')
                elif(len(x)-len(error_y)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(len(x)) + ' is not compatible with data of length ' + str(len(error_y)))
                else:
                    self.cov_y_plot = np.diag(self.error_y ** 2)
                    self.error_y = self.error_y[self.mask]
                    self.cov_y = np.diag(self.error_y**2)
                    self.inv_cov_y = np.abs(np.linalg.inv(self.cov_y))

            elif (error_x.any() != None and error_y.any() != None):
                if len(error_y.shape) > 1:
                    if error_y.shape[0] == 2 and isinstance(error_y.shape[1], int):
                        if (len(x) - len(error_y[0])) != 0 or (len(x) - len(error_y[1])) != 0:
                            raise SyntaxError('Length of data does not match. Data of length ' + str(
                                len(x)) + ' is not compatible with data of length ' + str(error_y.shape))
                        elif (len(error_x) - len(error_y[0])) != 0 or (len(error_x) - len(error_y[1])) != 0:
                            raise SyntaxError('Length of data does not match. Data of length ' + str(
                                len(error_x)) + ' is not compatible with data of length ' + str(error_y.shape))
                        else:
                            self.cov_x_plot = np.diag(self.error_x ** 2)
                            self.error_x = self.error_x[self.mask]
                            self.cov_x = np.diag(self.error_x**2)
                            self.cov_y_plot = np.array([np.diag(self.error_y[0]**2), np.diag(self.error_y[1]**2)])
                            self.error_y = self.error_y[:, self.mask]
                            self.cov_y = np.array([np.diag(self.error_y[0]**2), np.diag(self.error_y[1]**2)])
                            self.inv_cov_y = np.array(
                                [np.abs(np.linalg.inv(self.cov_y[0])), np.abs(np.linalg.inv(self.cov_y[1]))])
                    else:
                        raise SyntaxError(
                            'Wrong shape for error_y ' + str(error_y.shape) + ' must be of shape (2,n) or (n,)')
                elif (len(error_x)-len(error_y)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(len(error_x)) + ' is not compatible with data of length ' + str(len(error_y)))
                elif (len(x)-len(error_y)) != 0:
                    raise SyntaxError('Length of data does not match. Data of length ' + str(len(x)) + ' is not compatible with data of length ' + str(len(error_y)))
                else:
                    self.cov_x_plot = np.diag(self.error_x ** 2)
                    self.error_x = self.error_x[self.mask]
                    self.cov_x = np.diag(self.error_x**2)
                    self.cov_y_plot = np.diag(self.error_y**2)
                    self.error_y = self.error_y[self.mask]
                    self.cov_y = np.diag(self.error_y**2)
                    self.inv_cov_y = np.abs(np.linalg.inv(self.cov_y))

            else:
                raise NotImplementedError("No errors on x or y?")

        self.func = function
        self.jac = jac
        self.hess = hess
        self.x_plot = x
        self.y_plot = y
        self.x = x[self.mask]
        self.y = y[self.mask]

        self.method = method
        '''Check compatibility of P0'''
        if method != 'diff_evol':
            if P0.any() == None:
                self.params = np.ones(len(signature(function).parameters)-1)
                self.initparams = self.params
            else:
                if isinstance(P0, (np.ndarray)) != 1:
                    raise SyntaxError('P0 has to be a numpy array')
                elif (len(signature(function).parameters)-1) != len(P0):
                    raise SyntaxError('Length of P0 ' + str(len(P0)) + ' is not compatible with ' + str(len(signature(function).parameters)-1) + ' parameters in function')
                else:
                    self.initparams = P0
                    self.params = P0
        else:
            if P0.any() == None:
                self.params = np.ones((len(signature(function).parameters)-1,2))
                self.initparams = np.ones(len(signature(function).parameters))
            else:
                if isinstance(P0, (np.ndarray)) != 1:
                    raise SyntaxError('P0 has to be a numpy array')
                elif P0.shape != (len(signature(function).parameters)-1,2):
                    raise SyntaxError('Shape of P0 ' + str(P0.shape) + ' is not compatible with shape needed ' + str((len(signature(function).parameters)-1,2)))
                else:
                    self.initparams = P0[:,0]
                    self.params = P0

        self.df = len(x)-len(self.params)


        '''Fitting with Minimize'''
        if self.jac != None:
            if method != 'diff_evol':
                self.FitResult = minimize(self.Chi2, self.params,
                                          args=(self.func, self.jac, self.hess, self.x, self.y, self.cov_y, self.cov_x), jac=self.Chi2Deriv,
                                          method=method, **kwargs)
                self.params = self.FitResult.x
                self.Chi2 = self.FitResult.fun
                self.Cov = self.FitResult.hess_inv * 2
                self.Error = np.sqrt(np.diag(self.Cov))
                self.Pval = chi2.sf(self.Chi2, self.df)
                self.NormRes = self._NormRes(self.params, self.func, self.x, self.y, self.cov_y, self.cov_x)
            else:
                '''Fitting with Minimize and finding global minimum'''
                self.FitResult = differential_evolution(self.Chi2, self.params,
                                                args=(self.func, self.jac, self.hess, self.x, self.y, self.cov_y, self.cov_x), **kwargs)
                self.params = self.FitResult.x
                if self.hess != None:
                    self.Chi2 = self.FitResult.fun
                    self.Cov = np.linalg.inv(self.Chi2Deriv2(self.params,self.func,self.jac,self.hess,self.x,self.y,self.cov_y,self.cov_x))*2
                    self.Error = np.sqrt(np.abs(np.diag(self.Cov)))
                    self.Pval = chi2.sf(self.Chi2, self.df)
                    self.NormRes = self._NormRes(self.params, self.func, self.x, self.y, self.cov_y, self.cov_x)
                    if self.FitResult.success != 1:
                        print(self.FitResult.message + ' P-value at ' + str(self.Pval))
                else:
                    self.FitResult = minimize(self.Chi2, self.params,
                                              args=(self.func, self.jac, self.hess, self.x, self.y, self.cov_y, self.cov_x), jac=self.Chi2Deriv, hess=self.Chi2Deriv2)
                    self.params = self.FitResult.x
                    self.Chi2 = self.FitResult.fun
                    self.Cov = self.FitResult.hess_inv * 2
                    self.Error = np.sqrt(np.diag(self.Cov))
                    self.Pval = chi2.sf(self.Chi2, self.df)
                    self.NormRes = self._NormRes(self.params, self.func, self.x, self.y, self.cov_y, self.cov_x)
                    if self.FitResult.success != 1:
                        print(self.FitResult.message + ' P-value at ' + str(self.Pval))
        else:
            if method != 'diff_evol':
                self.FitResult = minimize(self.Chi2, self.params,
                                          args=(self.func, self.jac, self.hess, self.x, self.y, self.cov_y, self.cov_x),
                                          method=method, **kwargs)
                self.params = self.FitResult.x
                self.Chi2 = self.FitResult.fun
                self.Cov = self.FitResult.hess_inv * 2
                self.Error = np.sqrt(np.diag(self.Cov))
                self.Pval = chi2.sf(self.Chi2, self.df)
                self.NormRes = self._NormRes(self.params, self.func, self.x, self.y, self.cov_y, self.cov_x)
            else:
                '''Fitting with Minimize and finding global minimum'''
                params = differential_evolution(self.Chi2, self.params,
                                                args=(self.func, self.jac, self.hess, self.x, self.y, self.cov_y, self.cov_x), **kwargs).x
                self.FitResult = minimize(self.Chi2, params,
                                          args=(self.func, self.jac, self.hess, self.x, self.y, self.cov_y, self.cov_x))
                self.params = self.FitResult.x
                self.Chi2 = self.FitResult.fun
                self.Cov = self.FitResult.hess_inv * 2
                self.Error = np.sqrt(np.diag(self.Cov))
                self.Pval = chi2.sf(self.Chi2, self.df)
                self.NormRes = self._NormRes(self.params, self.func, self.x, self.y, self.cov_y, self.cov_x)
                if self.FitResult.success != 1:
                    print(self.FitResult.message + ' P-value at ' + str(self.Pval))


    def Plot(self, ms=5, capsize = 3, elinewidth=None, lw=None, xlabel='x', ylabel='y', figsize=(6,6), init=False, legend=True, stats=True, do_return=False, legend_loc = 'upper right', save=None):
        fig, axes = plt.subplots(2,1, sharex=True, figsize=figsize, gridspec_kw={'height_ratios': [2, 1]})

        if (self.cov_x_plot.any() == None and self.cov_y_plot.any() == None):
            axes[0].plot(self.x_plot,self.y_plot,'ko', ms=ms, capsize=capsize)
        elif (self.cov_x_plot.any() == None and self.cov_y_plot.any() != None):
            if len(self.cov_y_plot.shape) == 3:
                axes[0].errorbar(self.x_plot, self.y_plot, yerr=np.sqrt(np.array([np.diag(self.cov_y_plot[0]),np.diag(self.cov_y_plot[1])])),
                                 fmt='ko', ms=ms, capsize=capsize, elinewidth=elinewidth)
            else:
                axes[0].errorbar(self.x_plot,self.y_plot,yerr=np.sqrt(np.diag(self.cov_y_plot)),
                                 fmt='ko', ms=ms, capsize=capsize, elinewidth=elinewidth)
        elif (self.cov_x_plot.any() != None and self.cov_y_plot.any() == None):
            axes[0].errorbar(self.x_plot,self.y_plot,xerr=np.sqrt(np.diag(self.cov_x_plot)),
                             fmt='ko', ms=ms, capsize=capsize, elinewidth=elinewidth)
        elif (self.cov_x_plot.any() != None and self.cov_y_plot.any() != None):
            if len(self.cov_y_plot.shape) == 3:
                axes[0].errorbar(self.x_plot, self.y_plot, yerr=np.sqrt(np.array([np.diag(self.cov_y_plot[0]),np.diag(self.cov_y_plot[1])])),
                                 xerr=np.sqrt(np.diag(self.cov_x_plot)),
                                 fmt='ko', ms=ms, capsize=capsize, elinewidth=elinewidth)
            else:
                axes[0].errorbar(self.x_plot,self.y_plot,yerr=np.sqrt(np.diag(self.cov_y_plot)), xerr=np.sqrt(np.diag(self.cov_x_plot)),
                                 fmt='ko', ms=ms, capsize=capsize, elinewidth=elinewidth)
        else:
            print("Something is wrong")

        xp = np.linspace(np.min(self.x),np.max(self.x),1000,endpoint=True)
        axes[0].plot(xp,self.func(xp,*self.params),color='b', label = 'Fit', lw=lw)
        axes[0].set_ylabel(ylabel)
        if init and self.method != 'diff_evol':
            axes[0].plot(xp,self.func(xp,*self.initparams),color='r', label='Initial', lw=lw)
        if legend:
            if stats:
                handles, labels = axes[0].get_legend_handles_labels()
                patch = mpatches.Patch(linewidth=0,fill=False,label=f'P-val: {self.Pval:.2e}')
                handles.append(patch) 
                axes[0].legend(handles=handles,loc = legend_loc)
            else:
                axes[0].legend(loc = legend_loc)


        '''Residual plot'''
        axes[1].plot(self.x,self.NormRes, 'ko', ms=ms, lw=lw)
        axes[1].hlines(1,np.min(self.x),np.max(self.x), color = 'r', linestyle='dashed', lw=lw)
        axes[1].hlines(-1,np.min(self.x),np.max(self.x), color = 'r', linestyle='dashed', lw=lw)
        axes[1].set_ylabel('Norm. Res.')
        axes[1].set_xlabel(xlabel)

        plt.tight_layout()
        if save != None:
            plt.savefig(save)
        if do_return:
            return fig, axes
        else:
            plt.show()

    def derivative(self, params, func, x):
        eps = np.sqrt(np.finfo(float).eps) * (1.0 + x)
        return (func(x + eps, *params) - func(x - eps, *params)) / (2.0 * eps)

    def jacobian(self, params, func, x):
        matrix = np.empty((len(x), len(x)))
        for i in range(len(x)):
            for j in range(len(x)):
                matrix[i][j] = self.derivative(params, func, x[i]) * self.derivative(params, func, x[j])
        return matrix

    def jacobian_vec(self, params, func, x, jac=False, hess=False):
        if jac:
            vec = np.empty((len(x),len(params)))
            for i in range(len(x)):
                vec[i,:] = self.derivative(params, func, x[i]) * self.derivative(params, func, x[i])
        elif hess:
            vec = np.empty((len(x), len(params), len(params)))
            for i in range(len(x)):
                vec[i, :, :] = self.derivative(params, func, x[i]) * self.derivative(params, func, x[i])
        else:
            vec = np.empty(len(x))
            for i in range(len(x)):
                vec[i] = self.derivative(params, func, x[i]) * self.derivative(params, func, x[i])
        return vec

    def Chi2Deriv(self, params, func, param_jac, param_hess, x, y, cov_y, cov_x=np.array([None])):
        R = y - func(x, *params)
        F = param_jac(x, *params)
        Res_deriv = np.empty(len(params))
        if len(cov_y.shape) == 3:
            if cov_x.any() == None:
                for i in range(len(F)):
                    ResDerivLow = np.einsum('i,ij,j->i', F[i], self.inv_cov_y[0], R) + np.einsum('i,ij,j->i', R, self.inv_cov_y[0], F[i])
                    ResDerivHigh = np.einsum('i,ij,j->i', F[i], self.inv_cov_y[1], R) + np.einsum('i,ij,j->i', R, self.inv_cov_y[1], F[i])
                    ResDeriv = self.ResLowHigh(ResDerivLow, ResDerivHigh, R)
                    Res_deriv[i] = -np.sum(ResDeriv)
            else:
                if self.IsDiag(cov_y[0]) and self.IsDiag(cov_x) and self.IsDiag(cov_y[1]):
                    jac_f = self.jacobian_vec(params, func, x)
                    jac_F = self.jacobian_vec(params, param_jac, x, jac=True)
                    inv_cov_low = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y[0]))
                    inv_cov_high = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y[1]))
                    for i in range(len(F)):
                        grad_inv_cov = inv_cov_low * np.diagonal(cov_x) * jac_F[:, i] * inv_cov_low
                        ResDerivLow = R * grad_inv_cov * R + F[i] * inv_cov_low * R + R * inv_cov_low * F[i]

                        grad_inv_cov = inv_cov_high * np.diagonal(cov_x) * jac_F[:, i] * inv_cov_high
                        ResDerivHigh = R * grad_inv_cov * R + F[i] * inv_cov_high * R + R * inv_cov_high * F[i]

                        ResDeriv = self.ResLowHigh(ResDerivLow, ResDerivHigh, R)
                        Res_deriv[i] = -np.sum(ResDeriv)
                else:
                    jac_f = self.jacobian(params, func, x)
                    jac_F = self.jacobian(params, param_jac, x, jac=True)
                    inv_cov_low = np.abs(np.linalg.inv(np.einsum('ij,ij->ij', cov_x, jac_f) + cov_y[0]))
                    inv_cov_high = np.abs(np.linalg.inv(np.einsum('ij,ij->ij', cov_x, jac_f) + cov_y[1]))
                    for i in range(len(F)):
                        grad_inv_cov = np.einsum('il,lk,lk,kj->ij', inv_cov_low, cov_x, jac_F[:, :, i], inv_cov_low)
                        ResDerivLow = np.einsum('i,ij,j->i', R, grad_inv_cov, R) + np.einsum('i,ij,j->i', F[i], inv_cov_low,R) + \
                                 np.einsum('i,ij,j->i', R, inv_cov_low, F[i])

                        grad_inv_cov = np.einsum('il,lk,lk,kj->ij', inv_cov_high, cov_x, jac_F[:, :, i], inv_cov_high)
                        ResDerivHigh = np.einsum('i,ij,j->i', R, grad_inv_cov, R) + np.einsum('i,ij,j->i', F[i], inv_cov_high, R) + \
                                 np.einsum('i,ij,j->i', R, inv_cov_high, F[i])

                        ResDeriv = self.ResLowHigh(ResDerivLow, ResDerivHigh, R)
                        Res_deriv[i] = -np.sum(ResDeriv)
        else:
            if cov_x.any() == None:
                for i in range(len(F)):
                    Res_deriv[i] = np.nan_to_num(-np.sum(np.einsum('i,ij,j->i', F[i], self.inv_cov_y, R) + np.einsum('i,ij,j->i', R, self.inv_cov_y, F[i])),nan=1e200)
            else:
                if self.IsDiag(cov_y) and self.IsDiag(cov_x):
                    jac_f = self.jacobian_vec(params, func, x)
                    jac_F = self.jacobian_vec(params, param_jac, x, jac=True)
                    inv_cov = 1 / np.abs(np.diagonal(cov_x)*jac_f + np.diagonal(cov_y))
                    for i in range(len(F)):
                        grad_inv_cov = inv_cov * np.diagonal(cov_x) * jac_F[:,i] * inv_cov
                        Res_deriv[i] = -np.sum(R*grad_inv_cov*R + F[i]*inv_cov*R + R*inv_cov*F[i])
                else:
                    jac_f = self.jacobian(params, func, x)
                    jac_F = self.jacobian(params, param_jac, x, jac=True)
                    inv_cov = np.abs(np.linalg.inv(np.einsum('ij,ij->ij', cov_x, jac_f) + cov_y))
                    for i in range(len(F)):
                        grad_inv_cov = np.einsum('il,lk,lk,kj->ij', inv_cov, cov_x, jac_F[:,:,i], inv_cov)
                        Res_deriv[i] = -np.sum(np.einsum('i,ij,j->i', R, grad_inv_cov, R) + np.einsum('i,ij,j->i', F[i], inv_cov, R) + np.einsum(
                            'i,ij,j->i', R, inv_cov, F[i]))
        return Res_deriv

    def Chi2Deriv2(self, params, func, param_jac, param_hess, x, y, cov_y, cov_x=np.array([None])):
        R = y - func(x, *params)
        F = param_jac(x, *params)
        H = param_hess(x, *params)
        Res_deriv = np.empty((len(params), len(params)))
        if len(cov_y.shape) == 3:
            if cov_x.any() == None:
                for i in range(len(F)):
                    for j in range(len(F)):
                        ResDerivLow = np.einsum('i,ij,j->i', H[i, j], self.inv_cov_y[0], R) \
                                        - np.einsum('i,ij,j->i', F[i], self.inv_cov_y[0], F[j])
                        ResDerivHigh = np.einsum('i,ij,j->i', H[i, j], self.inv_cov_y[1], R) \
                                        - np.einsum('i,ij,j->i', F[i], self.inv_cov_y[1], F[j])

                        ResDeriv = self.ResLowHigh(ResDerivLow, ResDerivHigh, R)
                        Res_deriv[i,j] = -2*np.sum(ResDeriv)
            else:
                if self.IsDiag(cov_y[0]) and self.IsDiag(cov_x) and self.IsDiag(cov_y[1]):
                    jac_f = self.jacobian_vec(params, func, x)
                    jac_F = self.jacobian_vec(params, param_jac, x, jac=True)
                    jac_H = self.jacobian_vec(params, param_hess, x, hess=True)
                    inv_cov_low = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y[0]))
                    inv_cov_high = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y[1]))

                    grad_inv_cov_low = np.array([inv_cov_low * np.diagonal(cov_x) * jac_F[:, i] * inv_cov_low for i in range(len(F))])
                    grad_inv_cov_high = np.array([inv_cov_high * np.diagonal(cov_x) * jac_F[:, i] * inv_cov_high for i in range(len(F))])
                    for i in range(len(F)):
                        for j in range(len(F)):
                            gradsq_inv_cov = -grad_inv_cov_low[i] * np.diagonal(cov_x) * jac_F[:, j] * inv_cov_low \
                                             - inv_cov_low * np.diagonal(cov_x) * jac_F[:, j] * grad_inv_cov_low[i] \
                                             - inv_cov_low * np.diagonal(cov_x) * jac_H[:, i, j] * inv_cov_low
                            ResDerivLow = R * gradsq_inv_cov * R + \
                                          2 * (F[i] * grad_inv_cov_low[j] * R -
                                               F[j] * grad_inv_cov_low[i] * R +
                                               H[i, j] * inv_cov_low * R -
                                               F[i] * inv_cov_low * F[j])

                            gradsq_inv_cov = -grad_inv_cov_high[i] * np.diagonal(cov_x) * jac_F[:, j] * inv_cov_high \
                                             - inv_cov_high * np.diagonal(cov_x) * jac_F[:, j] * grad_inv_cov_high[i] \
                                             - inv_cov_high * np.diagonal(cov_x) * jac_H[:, i, j] * inv_cov_high
                            ResDerivHigh = R * gradsq_inv_cov * R + \
                                          2 * (F[i] * grad_inv_cov_high[j] * R -
                                               F[j] * grad_inv_cov_high[i] * R +
                                               H[i, j] * inv_cov_high * R -
                                               F[i] * inv_cov_high * F[j])

                            ResDeriv = self.ResLowHigh(ResDerivLow, ResDerivHigh, R)
                            Res_deriv[i, j] = -np.sum(ResDeriv)
                else:
                    jac_f = self.jacobian(params, func, x)
                    jac_F = self.jacobian(params, param_jac, x, jac=True)
                    jac_H = self.jacobian(params, param_hess, x, hess=True)

                    inv_cov_low = np.abs(np.linalg.inv(np.einsum('ij,ij->ij', cov_x, jac_f) + cov_y[0]))
                    inv_cov_high = np.abs(np.linalg.inv(np.einsum('ij,ij->ij', cov_x, jac_f) + cov_y[1]))

                    grad_inv_cov_low = np.array(
                        [[np.einsum('il,lk,lk,kj->ij', inv_cov_low, cov_x, jac_F[:, :, i], inv_cov_low)] for i in
                         range(len(F))])
                    grad_inv_cov_high = np.array(
                        [[np.einsum('il,lk,lk,kj->ij', inv_cov_high, cov_x, jac_F[:, :, i], inv_cov_high)] for i in
                         range(len(F))])

                    for i in range(len(F)):
                        for j in range(len(F)):
                            gradsq_inv_cov = -np.einsum('il,lk,lk,kj->ij', grad_inv_cov_low[i], cov_x, jac_F[:, :, j],
                                                        inv_cov_low) \
                                             - np.einsum('il,lk,lk,kj->ij', inv_cov_low, cov_x, jac_F[:, :, j],
                                                         grad_inv_cov_low[i]) \
                                             - np.einsum('il,lk,lk,kj->ij', inv_cov_low, cov_x, jac_H[:, :, i, j],
                                                         inv_cov_low)

                            ResDerivLow = np.einsum('i,ij,j->i', R, gradsq_inv_cov, R) + \
                                                      2 * (np.einsum('i,ij,j->i', F[i], grad_inv_cov_low[j], R) -
                                                           np.einsum('i,ij,j->i', F[j], grad_inv_cov_low[i], R) +
                                                           np.einsum('i,ij,j->i', H[i, j], inv_cov_low, R) -
                                                           np.einsum('i,ij,j->i', F[i], inv_cov_low, F[j]))

                            gradsq_inv_cov = -np.einsum('il,lk,lk,kj->ij', grad_inv_cov_high[i], cov_x, jac_F[:, :, j],
                                                        inv_cov_high) \
                                             - np.einsum('il,lk,lk,kj->ij', inv_cov_high, cov_x, jac_F[:, :, j],
                                                         grad_inv_cov_high[i]) \
                                             - np.einsum('il,lk,lk,kj->ij', inv_cov_high, cov_x, jac_H[:, :, i, j],
                                                         inv_cov_high)

                            ResDerivHigh = np.einsum('i,ij,j->i', R, gradsq_inv_cov, R) + \
                                          2 * (np.einsum('i,ij,j->i', F[i], grad_inv_cov_high[j], R) -
                                               np.einsum('i,ij,j->i', F[j], grad_inv_cov_high[i], R) +
                                               np.einsum('i,ij,j->i', H[i, j], inv_cov_high, R) -
                                               np.einsum('i,ij,j->i', F[i], inv_cov_high, F[j]))

                            ResDeriv = self.ResLowHigh(ResDerivLow, ResDerivHigh, R)
                            Res_deriv[i] = -np.sum(ResDeriv)
        else:
            if cov_x.any() == None:
                for i in range(len(F)):
                    for j in range(len(F)):
                        Res_deriv[i,j] = -2*np.sum(np.einsum('i,ij,j->i', H[i,j], self.inv_cov_y, R) -
                                               np.einsum('i,ij,j->i', F[i], self.inv_cov_y, F[j]))
            elif self.IsDiag(cov_y) and self.IsDiag(cov_x):
                    jac_f = self.jacobian_vec(params, func, x)
                    jac_F = self.jacobian_vec(params, param_jac, x, jac=True)
                    jac_H = self.jacobian_vec(params, param_hess, x, hess=True)
                    inv_cov = 1 / np.abs(np.diagonal(cov_x)*jac_f + np.diagonal(cov_y))
                    grad_inv_cov = np.array([inv_cov * np.diagonal(cov_x) * jac_F[:, i] * inv_cov for i in range(len(F))])
                    for i in range(len(F)):
                        for j in range(len(F)):
                            gradsq_inv_cov = -grad_inv_cov[i] * np.diagonal(cov_x) * jac_F[:, j] * inv_cov \
                                             -inv_cov * np.diagonal(cov_x) * jac_F[:, j] * grad_inv_cov[i]\
                                             -inv_cov * cov_x * jac_H[:, i, j] * inv_cov

                            Res_deriv[i,j] = -np.sum(R * gradsq_inv_cov * R +
                                                 2 * (F[i]*grad_inv_cov[j]*R -
                                                      F[j]*grad_inv_cov[i]*R +
                                                      H[i,j]*inv_cov*R -
                                                      F[i]*inv_cov*F[j]))
            else:
                jac_f = self.jacobian(params, func, x)
                jac_F = self.jacobian(params, param_jac, x, jac=True)
                jac_H = self.jacobian(params, param_hess, x, hess=True)

                inv_cov = np.abs(np.linalg.inv(np.einsum('ij,ij->ij', cov_x, jac_f) + cov_y))

                grad_inv_cov = np.array(
                    [[np.einsum('il,lk,lk,kj->ij', inv_cov, cov_x, jac_F[:, :, i], inv_cov)] for i in range(len(F))])

                for i in range(len(F)):
                    for j in range(len(F)):
                        gradsq_inv_cov = -np.einsum('il,lk,lk,kj->ij', grad_inv_cov[i], cov_x, jac_F[:, :, j],
                                                          inv_cov) \
                                               - np.einsum('il,lk,lk,kj->ij', inv_cov, cov_x, jac_F[:, :, j],
                                                           grad_inv_cov[i]) \
                                               - np.einsum('il,lk,lk,kj->ij', inv_cov, cov_x, jac_H[:, :, i, j],
                                                           inv_cov)

                        Res_deriv[i,j] = -np.sum(np.einsum('i,ij,j->i', R, gradsq_inv_cov, R) +
                                                 2 * (np.einsum('i,ij,j->i', F[i], grad_inv_cov[j], R) -
                                                      np.einsum('i,ij,j->i', F[j], grad_inv_cov[i], R) +
                                                      np.einsum('i,ij,j->i', H[i,j], inv_cov, R) -
                                                      np.einsum('i,ij,j->i', F[i], inv_cov, F[j])))
        return Res_deriv

    def Chi2(self, params, func, param_jac, param_hess, x, y, cov_y, cov_x=np.array([None])):
        Res = self.NormResSq(params, func, x, y, cov_y, cov_x)
        return np.sum(np.abs(Res))

    def NormResSq(self, params, func, x, y, cov_y, cov_x=np.array([None])):
        R = y - func(x, *params)
        if len(cov_y.shape) == 3:
            if cov_x.any() == None:
                ResLow = np.einsum('i,ij,j->i', R, self.inv_cov_y[0], R)
                ResHigh = np.einsum('i,ij,j->i', R, self.inv_cov_y[1], R)
                Res = self.ResLowHigh(ResLow, ResHigh, R)
            else:
                if self.IsDiag(cov_y[0]) and self.IsDiag(cov_x) and self.IsDiag(cov_y[1]):
                    jac_f = self.jacobian_vec(params, func, x)
                    inv_cov = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y[0]))
                    ResLow = R * inv_cov * R
                    inv_cov = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y[1]))
                    ResHigh = R * inv_cov * R
                    Res = self.ResLowHigh(ResLow, ResHigh, R)
                else:
                    jac_f = self.jacobian(params, func, x)
                    inv_cov = np.abs(np.linalg.inv(np.einsum('ij,jk->ik', cov_x, jac_f) + cov_y[0]))
                    ResLow = np.einsum('i,ij,j->i', R, inv_cov, R)
                    inv_cov = np.abs(np.linalg.inv(np.einsum('ij,jk->ik', cov_x, jac_f) + cov_y[1]))
                    ResHigh = np.einsum('i,ij,j->i', R, inv_cov, R)
                    Res = self.ResLowHigh(ResLow, ResHigh, R)
        else:
            if cov_x.any() == None:
                Res = np.einsum('i,ij,j->i', R, self.inv_cov_y, R)
            else:
                if self.IsDiag(cov_y) and self.IsDiag(cov_x):
                    jac_f = self.jacobian_vec(params, func, x)
                    inv_cov = 1 / np.abs(np.diagonal(cov_x) * jac_f + np.diagonal(cov_y))
                    Res = R * inv_cov * R
                else:
                    jac_f = self.jacobian(params, func, x)
                    inv_cov = np.abs(np.linalg.inv(np.einsum('ij,jk->ik', cov_x, jac_f) + cov_y))
                    Res = np.einsum('i,ij,j->i', R, inv_cov, R)
        return Res

    def _NormRes(self, params, func, x, y, cov_y, cov_x=np.array([None])):
        R = y - func(x, *params)
        Res = self.NormResSq(params, func, x, y, cov_y, cov_x)
        return np.array([np.sign(i)*np.abs(j) for i,j in zip(R,Res)])

    def ResLowHigh(self, ResLow, ResHigh, R):
        Res = np.empty(len(R))
        for i in range(len(Res)):
            if np.sign(R[i]) == 1:
                Res[i] = ResHigh[i]
            else:
                Res[i] = ResLow[i]
        return Res

    def IsDiag(self, M):
        i, j = M.shape
        assert i == j
        test = M.reshape(-1)[:-1].reshape(i - 1, j + 1)
        return ~np.any(test[:, 1:])

# %%
