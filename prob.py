from scipy import stats as st
import pyDOE
import math
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
from numpy import linspace
from numpy import var

# TODO figure out how to avoid negative numbers when sampling from some distrubtions/ results get messed up in some cases

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

        #print(g_norm_var)
        #print(g_norm_mean)

        g_lognorm_mu = 1#math.exp(g_norm_mean+(math.pow(g_norm_var,2)/2))
        g_lognorm_sigma = .5#math.exp(2*g_norm_mean+math.pow(g_norm_var,2))*math.exp(math.pow(g_norm_var,2)-1)

        g_uni_b = 2*math.pow(g_norm_var,2)
        g_uni_a = g_norm_mean - math.pow(g_norm_var,2)

        return [[g_norm_mean, g_norm_var],[g_lognorm_mu,g_lognorm_sigma], [g_uni_a,g_uni_b]]


    def make_pdf_hist(self) -> list:


        hist = {}

        bins, step = linspace(self.values[0],self.values[len(self.values)-1], self.num_bins, retstep=True)
        kde = stats.gaussian_kde(self.values)


        return [bins, kde]

    def plot_info(self) -> None:


        fit_params = self.cdfs[self.index][1]

        fig, ax = plt.subplots(1, 2, figsize=(12, 6))

        # histogram
        # hist_plot = list(sorted(self.hist.items(), key=lambda x: x[0]))
        # x = [i[0] for i in hist_plot]
        # y = [i[1] for i in hist_plot]
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

    def make_cdfs(self) -> (list,list):

        width = 1/len(self.values)
        y = []
        for i in range (len(self.values)):
            y.append(0+(width*i))

        cdf_list = [[stats.norm.cdf], [stats.lognorm.cdf], [stats.uniform.cdf]]

        for i in range (len(cdf_list)):
            param = optimize.curve_fit(cdf_list[i][0], self.values, y, p0=self.init_guess[i])[0]

            cdf_list[i].append(param)

        return cdf_list, y


    def ks_cdf(self) -> int:

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

        self.count += 1
        string = str(self.dist_types[self.index] + '(')

        for i in range (len(self.cdfs[self.index][1])-1):
            string += str(self.cdfs[self.index][1][i])
            string += ', '
        string += str(self.cdfs[self.index][1][i+1])
        string += ')'
        return string


def set_hyper_cube(model_para, Var):

    v_iter = int(model_para[0])
    u_iter = int(model_para[1])
    bin_num = int(model_para[2])

    if Var.type == 'V':
        hype_sample = v_iter//bin_num + 10
    else:
        hype_sample = u_iter//bin_num + 10

    lhs = pyDOE.lhs(bin_num, samples=hype_sample)
    lhs = lhs.ravel()
    Var.lhs = lhs


def sample_dist(name, Var, i, type):

    if Var.type == type or type == 'both':
        if Var.dist == 'Normal':
            mean = Var.param[0]
            std = Var.param[1]
            point = Var.lhs[i]
            sample = st.norm(loc=mean, scale=std).ppf(point)

        elif Var.dist == 'Uniform':
            a = Var.param[0]
            b = Var.param[1]
            point = Var.lhs[i]
            sample = st.uniform(loc=a, scale=b).ppf(point)

        elif Var.dist == 'Triangle':
            a = Var.param[0]
            b = Var.param[1]
            c = (Var.param[2]-a)/b
            point = Var.lhs[i]
            sample = st.triang(c, loc=a, scale=b).ppf(point)

        elif Var.dist == 'Log-Normal':
            mu_log = Var.param[0]
            sigma_log = Var.param[1]
            mu, sigma = lognorm_to_norm(mu_log,sigma_log)
            point = Var.lhs[i]
            sample = st.lognorm(sigma, scale=mu).ppf(point)

        elif Var.dist == 'Log-Uniform':
            a = Var.param[0]
            b = Var.param[1]
            point  = Var.lhs[i]
            sample = loguniform(loc=a, scale=b).ppf(point)

        else:
            print(Var.dist)
            print('There is a unknown distribution in', '\''+name + '\'')
            exit(0)


        Var.value = float(sample)



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

