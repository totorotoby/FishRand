import scipy.stats as stats
import matplotlib.pyplot as plt
from numpy import linspace
import numpy
import math
import scipy
from scipy import optimize


def my_log_normal_cdf(x, M_x, Sig_x):
    to_run = (math.log(x) - M_x) / Sig_x
    F_x = stats.norm.cdf(to_run)

    print('mine: ', F_x)


def fit_test():

    samples = lognormal_samples(0,.25)
    samples = sorted(samples)
    y = linspace(0, 1,num=len(samples))
    s, scale = lognorm_to_input(2, 3)
    m_x = math.log(scale)
    sig_x = s

    param = optimize.curve_fit(my_log_normal_cdf, samples, y, p0=[m_x,sig_x])[0]
    print(param)
    scipy_param = math.exp(param[0])
    print(scipy_param)
    total_width = max(samples) - min(samples)
    interval = total_width / 1000

    x1 = []
    for i in range(1000):
        x1.append(interval * i)

    plt.plot(samples,y, 'ro')
    plt.plot(x1, stats.lognorm.cdf(x1,s=param[1], scale=scipy_param))
    plt.show()


def my_log_normal_cdf(x, M_y, Sig_y):
    to_run = (numpy.log(x) - M_y) / Sig_y
    F_x = stats.norm.cdf(to_run)

    return F_x



def lognorm_to_input(M_y,Sig_y):
    print(M_y,Sig_y)
    m_x = (2 * math.log(M_y)) - (.5) * (math.log(math.pow(Sig_y, 2) + math.pow(M_y, 2)))

    scale = math.exp(m_x)

    sigma2 = -2 * math.log(M_y) + math.log(math.pow(Sig_y, 2) + math.pow(M_y, 2))
    s = math.sqrt(sigma2)

    return s, scale


def lognormal_test_of_ks_test():

    samples, s, scale = lognormal_samples(1, .25)

    my_args = [s, 0, scale]
    ks = stats.kstest(samples, 'lognorm', args=my_args)[0]


    #ks = stats.kstest(samples, 'lognorm', args=[s,scale])[0]
    print('ks: ', ks)

def lognormal_samples(M_y, Sig_y):

    m_x = (2*math.log(M_y)) - (.5)*(math.log(math.pow(Sig_y,2) + math.pow(M_y,2)))

    scale = math.exp(m_x)

    sigma2 = -2 * math.log(M_y) + math.log(math.pow(Sig_y,2) + math.pow(M_y,2))
    s = math.sqrt(sigma2)


    result = stats.lognorm(s, scale=scale).rvs(size=10000)

    return result, s, scale

