from scipy import stats as st
import pyDOE
import math
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
from numpy import linspace
from numpy import var

class Var:

    def __init__(self, type, dist, param):

        self.type = type   # 0 is uncertain , 1 is variable, 2 is both, 3 is point estimate
        self.dist = dist   # distribution name
        self.param = param
        self.values = None
        self.lhs = None

    def __str__(self):

        to_print = '\ntype | ' + self.type + '\ndistribution | ' + self.dist +'\nparameters | ' + str(self.param)
        return to_print


    def take_samples(self):

        if self.dist == 'Normal':
            mean = self.param[0]
            std = self.param[1]
            self.values = st.norm(loc=mean, scale=std).ppf(self.lhs)
        elif self.dist == 'Uniform':
            a = self.param[0]
            b = self.param[1]
            self.values = st.uniform(loc=a, scale=b).ppf(self.lhs)
        elif self.dist == 'Triangle':
            a = self.param[0]
            b = self.param[1]
            c = (self.param[2] - a) / b
            self.values = st.triang(c, loc=a, scale=b).ppf(self.lhs)
        elif self.dist == 'Log-Normal':
            mu_log = self.param[0]
            sigma_log = self.param[1]
            mu, sigma = lognorm_to_norm(mu_log,sigma_log)
            self.values = st.lognorm(sigma, scale=mu).ppf(self.lhs)
        elif self.dist == 'Log-Uniform':
            a = self.param[0]
            b = self.param[1]
            self.values = loguniform(loc=a, scale=b).ppf(self.lhs)
        elif self.dist == 'Beta':
            alpha = self.param[0]
            beta =  self.param[1]
            self.values = st.beta(alpha,beta).ppf(self.lhs)
        elif self.dist == 'Weibull':
            lamb = self.param[0]
            k = self.param[1]
            self.values = st.weibull_min(c=k, scale = lamb).ppf(self.lhs)
        else:
            print(self.dist)
            print('There is a unknown distribution called ', '\''+ self.dist + '\'')
            exit(0)


class loguniform:

    def __init__(self, loc=-1, scale=0, base=math.e):
        self.loc = loc
        self.scale = scale
        self.base = base

    def ppf(self,q):

        uniform = st.uniform(loc=self.loc, scale=self.scale)
        return math.pow(self.base, uniform.ppf(q))




class ResultDist:

    dist_types = ['norm', 'lognorm', 'uniform']

    def __init__(self, values, chemical, animal):

        self.count = 0
        self.chem = chemical
        self.animal = animal
        self.values = values
        self.values.sort()
        self.init_guess = self.inital_guess()
        self.num_bins = len(self.values)//10
        self.hist = self.make_pdf_hist()
        self.cdfs, self.y = self.make_cdfs()
        self.index = self.ks_cdf()
        self.best_para = self.bestparam()
        #self.plot_info()



    def inital_guess(self):

        g_norm_mean = 0
        for value in self.values:
            g_norm_mean += value

        g_norm_mean = g_norm_mean / len(self.values)
        g_norm_var = var(self.values)
        g_norm_std = math.sqrt(g_norm_var)

        g_lognorm_mu = math.exp(g_norm_mean+(math.pow( g_norm_std,2)/2))
        g_lognorm_sigma = math.exp(2*g_norm_mean+math.pow(g_norm_std,2))*math.exp(math.pow(g_norm_std,2)-1)

        g_uni_b = 2*math.pow( g_norm_std,2)
        g_uni_a = g_norm_mean - math.pow( g_norm_std,2)

        return [[g_norm_mean, g_norm_std],[g_lognorm_mu,g_lognorm_sigma], [g_uni_a,g_uni_b]]


    def make_pdf_hist(self):


        hist = {}

        bins, step = linspace(self.values[0],self.values[len(self.values)-1], self.num_bins, retstep=True)
        kde = stats.gaussian_kde(self.values)


        return [bins, kde]

    def plot_info(self):


        fit_params = self.cdfs[self.index][1]

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))

        pdf_x_label = '(ng/g) of ' + self.chem + ' in ' + self.animal
        pdf_y_label = 'P(x)'
        cdf_y_label = 'P(X < x)'
        ax[0].title.set_text('Probability Density Function')
        ax[1].title.set_text('Cumulative Density Function')
        ax[0].set_xlabel(pdf_x_label, size='large')
        ax[0].set_ylabel(pdf_y_label, size='large')
        ax[1].set_xlabel(pdf_x_label, size='large')
        ax[1].set_ylabel(cdf_y_label, size='large')
        x = self.hist[0]
        xmin = min(self.hist[0])
        xmax = max(self.hist[0])

        if self.dist_types[self.index] == 'norm':
            mean = fit_params[0]
            std = fit_params[1]
            totalwidth = 2*(mean+(3*std))
            x1 = self.make_x1(totalwidth)

            ax[0].set_xlim(xmin, xmax)
            ax[0].plot(x, self.hist[1](x), 'ro', color='r')
            ax[0].plot(x1, stats.norm.pdf(x1, loc=mean, scale=std), linewidth=2.0)
            ax[1].set_xlim(self.values[0], self.values[len(self.values)-1])
            ax[1].plot(self.values,self.y, 'ro', color='r')
            ax[1].plot(x1, stats.norm.cdf(x1, loc=mean, scale=std), color='g')
            plt.show()

        if self.dist_types[self.index] == 'uniform':

            a = fit_params[0]
            b = fit_params[1]
            totalwidth = 2*(a+b)
            x1 = self.make_x1(totalwidth)

            ax[0].set_xlim(xmin, xmax)
            ax[0].plot(x, self.hist[1](x), 'ro', color='r')
            ax[0].plot(x1, stats.uniform.pdf(x1, loc=a, scale=b), linewidth=2.0)
            ax[1].set_xlim(self.values[0], self.values[len(self.values) - 1])
            ax[1].plot(self.values, self.y, 'ro', color='r')
            ax[1].plot(x1, stats.uniform.cdf(x1, loc=a, scale=b), color='g')
            plt.show()

        if self.dist_types[self.index] == 'lognorm':
            mu = fit_params[0]
            sig2 = fit_params[1]
            totalwidth = 2*(mu+(3*sig2))
            x1 = self.make_x1(totalwidth)

            ax[0].set_xlim(xmin, xmax)
            ax[0].plot(x, self.hist[1](x), 'ro', color='r')
            ax[0].plot(x1, stats.lognorm.pdf(x1, 1 ,loc=mu, scale=sig2), linewidth=2.0)
            ax[1].set_xlim(self.values[0], self.values[len(self.values) - 1])
            ax[1].plot(self.values, self.y, 'ro', color='r')
            ax[1].plot(x1, stats.lognorm.cdf(x1, 1, loc=mu, scale=sig2), color='g')
            plt.show()



    def make_x1(self,total_width):
        interval = total_width / 1000

        x1 = []
        for i in range(1000):
            x1.append(interval * i)

        return x1

    def make_cdfs(self):

        width = 1/len(self.values)
        y = []
        for i in range (len(self.values)):
            y.append(0+(width*i))

        cdf_list = [[stats.norm.cdf], [stats.lognorm.cdf], [stats.uniform.cdf]]

        for i in range (len(cdf_list)):
            param = optimize.curve_fit(cdf_list[i][0], self.values, y, p0=self.init_guess[i])[0]

            cdf_list[i].append(param)



        return cdf_list, y


    def ks_cdf(self):

        ks_list = []

        for i in range (len(self.dist_types)):
            ks = stats.kstest(self.values, self.dist_types[i], args=self.cdfs[i][1])[0]
            ks_list.append(ks)

        min = 2
        index = 0
        for i in range (len(ks_list)):
            if ks_list[i] < min:
                min = ks_list[i]
                index = i

        return index

    def bestparam(self):

        params = []
        self.count += 1
        string = str(self.animal) + ' ' + str(self.chem) + ' '
        string += str(self.dist_types[self.index] + '(')

        for i in range (len(self.cdfs[self.index][1])-1):
            params.append(self.cdfs[self.index][1][i])
            string += str(self.cdfs[self.index][1][i])
            string += ', '
        params.append(self.cdfs[self.index][1][i+1])
        string += str(self.cdfs[self.index][1][i+1])
        string += ')'
        return [string, params]


def set_hyper_samp_cube(model_para, Var):

    v_iter = int(model_para[0])
    u_iter = int(model_para[1])
    bin_num = int(model_para[2])

    if Var.type == 'V':
        hype_sample = v_iter*u_iter//bin_num + 2
    else:
        hype_sample = u_iter//bin_num + 2

    lhs = pyDOE.lhs(bin_num, samples=hype_sample)
    lhs = lhs.ravel()
    Var.lhs = lhs
    Var.take_samples()


def lognorm_to_norm(mu_log,sigma_log):

    mu = 2*math.log(mu_log) - (1/2)*math.log((sigma_log*sigma_log) + (mu_log*mu_log))
    sigma = math.sqrt(-2*(math.log(mu_log))+math.log((sigma_log*sigma_log) + (mu_log*mu_log)))

    return mu, sigma


def make_result_dist(dicts):

    result_dict = {}

    for region, values in dicts[0].items():
        result_dict[region] = {}
        for animal, values1 in values.items():
            result_dict[region][animal] = {}
            for chemical in values1.keys():
                values = get_values(region,chemical, animal, dicts)
                new_dist = ResultDist(values, chemical,animal)
                result_dict[region][animal][chemical] = new_dist


    return result_dict


def get_values(region,chem, animal, dicts):

    dist_list = []

    for dict in dicts:
        dist_list.append(dict[region][animal][chem])

    return dist_list

