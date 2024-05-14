"""
functions.py
============

This module contains all the functions necessary to reproduce the results reported
on the paper entitled "Ratchet current and scaling properties in a nontwist mapping",
by Matheus Rolim Sales, Daniel Borin, Leonardo Costa de Souza, José Danilo Szezech Jr.,
Ricardo Luiz Viana, Iberê Luiz Caldas, e Edson Denis Leonel.

Functions
---------

1. map(theta0: float, I0: float, eps: float) -> ndarray:
    Calculates one iteration of the mapping given the initial condition (`theta0`, `I0`).
    
2. time_series(theta0: float, I0: float, eps: float, N: int) -> ndarray:
    Calculates the time series of the mapping given an initial condition (`theta0`, `I0`).
    
3. av_sqr_action(theta0: ndarray, I0: ndarray, eps: float, N: int) -> ndarray:
    Calculates the square root of the averaged action, I_rms, as a function of time.
    
4. lyapunov(theta0: float, I0: float, eps: float, N: int) -> ndarray:
    Calculates the largest Lyapunov exponent.
    
5. recurrence_times(theta0: float, I0: float, eps: float, I_rec: int, N: int) -> ndarray:
    Calculates the recurrence times for the chosen region defined by `I_rec`.
    
6. cumul_rtd(rtd: ndarray, t_min: int = 1e0, t_max: int = 1e6, num_t: int = 250) -> Tuple[ndarray, ndarray]:
    Calculates the cumulative distribution of recurrence times given the recurrence time distribution.
    
7. return_moments(y: ndarray, ms: ndarray) -> ndarray:
    Calculates the higher moments of the distribution `y`.
    
8. escape_time(theta0: float, I0: float, eps: float, I_esc: float, Nmax: int) -> ndarray:
    Calculates the escape time of a given initial condition for the survival region defined by `I_esc`.
    
9. escape_basin(theta0: float, I0: float, eps: float, I_esc: float, Nmax: int) -> ndarray:
    Calculates the escape basin of a given initial condition for the survival region defined by `I_esc`.
    
10. escape_time_and_basin(theta0: float, I0: float, eps: float, I_esc: float, Nmax: int) -> ndarray:
    Calculates the escape time and basin of a given initial condition for the survival region defined by `I_esc`.
    
11. survival_prob(escape_times: ndarray, N: int) -> ndarray:
    Calculates the survival probability given the escape times.
    
12. fixed_points(I: ndarray, m: int) -> float:
    Equation to be solved to return the fixed points.
    
13. plot_params(fontsize: int = 20, legend_fontsize: int = 14, axes_linewidth: float = 1.3) -> Tuple[int, int, float]:
    Update the parameters of the plot.

Author: Matheus Rolim Sales
Email: matheusrolim95@gmail.com
Last updated: 08/05/2024
"""

from params import *
import numpy as np # NumPy module
from numba import vectorize, njit # Numba module to create fast functions
import matplotlib.pyplot as plt
import matplotlib as mpl

@njit
def map(theta0: np.float64, I0: np.float64, eps: np.float64) -> np.ndarray:
    """
    Calculates one iteration of the mapping given the initial condition (`theta0`, `I0`).

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.

    Return
    ------
    ndarray
        The next value for the coordinates of the map.
    """
    theta = theta0
    I = I0

    I = I + eps * np.sin(2*np.pi*theta)
    q_factor = q1 + q2 * I**2 + q3 * I**3
    v = v1 + v2 * np.tanh(v3 * I + v4)
    E = e1 * I + e2 * np.sqrt(np.abs(I)) + e3
    theta = (theta + mu * v * (M / q_factor - L) + rho * E / np.sqrt(np.abs(I))) % 1.0

    return np.array([theta, I])

@njit
def time_series(theta0: np.float64, I0: np.float64, eps: np.float64, N: np.int32) -> np.ndarray:
    """
    Calculates the time series of the mapping given an initial condtion (`theta0`, `I0`).

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.
    N : int
        The length of the times series.

    Return
    ------
    ndarray
        The time series of the variables theta and I, with u[:, 0] = theta and u[:, 1] = I.
    """
    u = np.zeros((N + 1, 2))

    theta = theta0
    I = I0

    u[0, 0] = theta0
    u[0, 1] = I0
    for j in range(1, N):
        
        theta, I = map(theta, I, eps)

        u[j, 0] = theta
        u[j, 1] = I

    return u

@njit
def av_sqr_action(theta0: np.ndarray, I0: np.ndarray, eps: np.float64, N: np.int32) -> np.ndarray:
    """
    Calculates the square root of the averaged action, I_rms, as a function of time.

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    theta0 : np.ndarray
        The initial angle of the ensemble of initial conditions.
    I0 : np.ndarray
        The initial action of the ensemble of initial conditions.
    eps : float
        The perturbation.
    N : int
        The length of the times series.

    Return
    ------
    ndarray
        The I_rms as a function of time.
    """
    n_part = theta0.shape[0]
    
    theta = theta0.copy()
    I = I0.copy()
    
    S = np.zeros(n_part, dtype=np.float64)
    Irms = np.zeros(N, dtype=np.float64)
    
    for i in range(1, N + 1):
        theta, I = map(theta, I, eps)

        S += I ** 2
        Irms[i - 1] = np.sqrt((sum(S)/i)/n_part)
    
    return Irms

@njit
def lyapunov(theta0: np.float64, I0: np.float64, eps: np.float64, N: np.int32):
    """
    Calculates the largest Lyapunov exponent.

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.
    N : int
        The length of the times series.

    Return
    ------
    ndarray
        The largest Lyapunov exponent and the final value of theta and I.
    """
    # Original IC
    theta = theta0
    I = I0
    # Perturbed IC
    theta_p = (theta + 1e-8)
    I_p = I

    lypnv = 0
    for i in range(N):
        # Initial distance between the original and perturbed IC
        delta0 = np.sqrt((theta - theta_p)**2 + (I - I_p)**2)
        # Iterate the original IC
        theta, I = map(theta, I, eps)
        # Iterate the perturbed IC
        theta_p, I_p = map(theta_p, I_p, eps)
        # Final distance between the original and perturbed IC
        delta1 = np.sqrt((theta - theta_p)**2 + (I - I_p)**2)
        #
        lypnv += np.log(delta1/delta0)
        # Rescales the perturbed IC
        theta_p = (theta + delta0*(theta_p - theta)/delta1) % 1.0
        I_p = I + delta0*(I_p - I)/delta1

    return np.array([lypnv/N, theta, I])

@njit
def recurrence_times(theta0: np.float64, I0: np.float64, eps: np.float64, I_rec:np.int32, N: np.int32) -> np.ndarray:
    """
    Calculates the recurrence times for the chosen region defined by I_rec.

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.
    I_rec : int
        The recurrence region. If I_rec = 1 (-1), calculates the recurrence times
        for the upper (lower) region.
    N : int
        The length of the times series.

    Return
    ------
    ndarray
        The recurrence times.
    """
    theta, I = theta0, I0
    rec_times = []
    rt = 0
    has_entered = False
    for i in range(N):
        theta, I = map(theta, I, eps)
        # Check if the IC has entered the recurrence region
        if I*I_rec > 0: # Check if I and I_old have different sign. If so, the IC has crossed the line I = 0
            has_entered = True
            rt += 1
        elif has_entered:
            has_entered = False
            rec_times.append(rt)
            rt = 0
    
    return np.array(rec_times)

@njit
def cumul_rtd(rtd, t_min=1e0, t_max=1e6, num_t = 250):
    """
    Calculates the cumulative distribution of recurrence times given the recurrence time distribution.

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    rtd : np.ndarray
        The recurrence times
    t_min : int, optional
        The minimum time for the cumulative distribution. Default = 10^0
    t_min : int, optional
        The minimum time for the cumulative distribution. Default = 10^6
    num_t : int, optional
        The size of the sample.

    Return
    ------
    ndarray
        The cumulative distribution of recurrence times in log scale.
    """
    rtd_sorted = np.sort(rtd)[::-1]
    N_t = np.zeros(num_t)
    t = np.logspace(np.log10(t_min), np.log10(t_max), num_t)
    for j in range(num_t):
        rec_t = rtd_sorted >= t[j]
        aux = np.where(rec_t == False)[0]
        if len(aux) == 0:
            aux = -1
        else:
            aux = aux[0]
        rec_t = rec_t[:aux]
        N_t[j] = len(rec_t)
    
    Q = N_t/len(rtd)

    return t, Q

def return_moments(y: np.ndarray, ms: np.ndarray) -> np.ndarray:
    """
    Calculates the higher moments of the distribution y.

    ----------
    y : np.ndarray
        The probability distribution.
    ms : np.ndarray
        An array with the order of the moments.
    
    Return
    ------
    ndarray
        The higher moments of y normalized to y^m.
    """
    y_mean = np.mean(y)
    ym_mean = np.zeros_like(ms, dtype=np.float128)
    for m in ms:
        ym_mean[m] = np.mean(y**m)/(y_mean**m)
    return ym_mean

@vectorize(["i8(f8, f8, f8, f8, i8)"],
           nopython=True,
           target="parallel")
def escape_time(theta0: np.float64, I0: np.float64, eps: np.float64, I_esc: np.float64, Nmax: np.int32) -> np.ndarray:
    """
    Calculates the escape time of a given initial condition for the survival region defined by `I_esc`.

    This function uses Numba's `vectorize` decorator for performance optimization and vectorization of the function.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.
    I_esc : float
        The limits of the survival region.
    Nmax : int
        The maximum number of iterations.
    
    Return
    ------
    int
        The time it takes for the given initial condtion to escape the survival region.
    """
    theta = theta0
    I = I0
    for i in range(1, Nmax + 1):
        theta, I = map(theta, I, eps)
        if I < -I_esc:
            break
        elif I > I_esc:
            break
    
    return i

@vectorize(["i8(f8, f8, f8, f8, i8)"],
           nopython=True,
           target="parallel")
def escape_basin(theta0: np.float64, I0: np.float64, eps: np.float64, I_esc: np.float64, Nmax: np.int32) -> np.ndarray:
    """
    Calculates the escape basin of a given initial condition for the survival region defined by `I_esc`.

    This function uses Numba's `vectorize` decorator for performance optimization and vectorization of the function.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.
    I_esc : float
        The limits of the survival region.
    Nmax : int
        The maximum number of iterations.
    
    Return
    ------
    int
        The exit through which the given initial condition escapes the survival region.
    """
    theta = theta0
    I = I0
    esc = 0
    for i in range(1, Nmax + 1):
        theta, I = map(theta, I, eps)
        if I < -I_esc:
            esc = -1
            break
        elif I > I_esc:
            esc = 1
            break
    
    return esc

@njit
def escape_time_and_basin(theta0: np.float64, I0: np.float64, eps: np.float64, I_esc: np.float64, Nmax: np.int32) -> np.ndarray:
    """
    Calculates the escape time and basin of a given initial condition for the survival region defined by `I_esc`.

    This function uses Numba's `njit` decorator for performance optimization.

    Parameters
    ----------
    theta0 : float
        The initial angle.
    I0 : float
        The initial action.
    eps : float
        The perturbation.
    I_esc : float
        The limits of the survival region.
    Nmax : int
        The maximum number of iterations.
    
    Return
    ------
    ndarray
        The time and exit through which the given initial condition escapes the survival region.
    """
    theta = theta0
    I = I0
    esc = 0
    for i in range(1, Nmax + 1):
        I = I + eps*np.sin(2*np.pi*theta)
        q_factor = q1 + q2*I**2 + q3*I**3
        v = v1 + v2*np.tanh(v3*I + v4)
        E = e1*I + e2*np.sqrt(np.abs(I)) + e3
        theta = (theta + mu*v*(M/q_factor - L) + rho*E/np.sqrt(np.abs(I))) % 1.0
        if I < -I_esc:
            esc = -1
            break
        elif I > I_esc:
            esc = 1
            break
    
    return np.array([i, esc])

@njit
def survival_prob(escape_times: np.ndarray, N: np.int32):
    """
    Calculates the survival probability given the escape times.

    Parameters
    ----------
    escape_times : np.ndarray
        The escape times.
    N : int
        The maximum number of iterations when evaluating the escape times.

    Return
    ------
    ndarray
        The survival probability as a function of time.
    """
    sp = np.zeros(N)
    n_ic = len(escape_times)

    escape_times.sort()

    for i in range(N):
        # Using binary search to find the first index where escape_times[j] > (i + 1)
        idx = np.searchsorted(escape_times, i + 1, side='right')
        sp[i] = n_ic - idx

    return sp / n_ic

def fixed_points(I: np.ndarray, m: np.int32):
    """
    Equation to be solved to return the fixed points

    Parameters
    ----------
    I : float
        The value of the action.
    m : int
        An integer.

    Return
    ------
    float
        Return the function F(I) = 0.
    """
    q_factor = q1 + q2 * I**2 + q3 * I**3
    v = v1 + v2 * np.tanh(v3 * I + v4)
    E = e1 * I + e2 * np.sqrt(abs(I)) + e3

    return mu * v * (M / q_factor - L) + rho * E / np.sqrt(abs(I)) - m

def plot_params(fontsize: np.int32 = 20, legend_fontsize: np.int32 = 14, axes_linewidth: np.float64 = 1.3) -> np.ndarray:
    """
    Update the parameters of the plot.

    Returns
    -------
    np.ndarray
        The parameters passed to the function.
    """
    tick_labelsize = fontsize - 3
    plt.clf()
    plt.rc('font', size=fontsize)
    plt.rc('xtick', labelsize=tick_labelsize)
    plt.rc('ytick', labelsize=tick_labelsize)
    plt.rc('legend', fontsize=legend_fontsize)
    font = {'family' : 'stix'}
    plt.rc('font', **font)
    plt.rcParams["mathtext.fontset"] = "stix"
    mpl.rcParams['axes.linewidth'] = axes_linewidth #set the value globally

    return fontsize, legend_fontsize, axes_linewidth