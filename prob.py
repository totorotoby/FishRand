from scipy import stats as st
import pyDOE
import math
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
from numpy import linspace
import numpy

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
            b = self.param[1] - self.param[0]
            c = (self.param[2] - a) / b
            self.values = st.triang(c, loc=a, scale=b).ppf(self.lhs)
        elif self.dist == 'Log-Normal':
            m_y = self.param[0]
            sig_y = self.param[1]
            s, scale = lognorm_to_scipyinput(m_y,sig_y)
            self.values = st.lognorm(s=s, scale=scale).ppf(self.lhs)
        elif self.dist == 'Log-Uniform':
            a = self.param[0]
            b = self.param[1]
            self.values = loguniform(loc=a, scale=b).ppf(self.lhs)
        elif self.dist == 'Beta':
            alpha = self.param[0]
            beta =  self.param[1]
            self.values = st.beta(alpha,beta).ppf(self.lhs)
            self.values = sorted(self.values)
            totalwidth = max
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
        #self.hist = self.make_pdf_hist()
        self.cdfs, self.y = self.make_cdfs()
        self.index = self.ks_cdf()
        #self.show_all_cdf()
        self.best_para = self.bestparam()
        #self.plot_info()

    def inital_guess(self):
        # guess for normal
        M_y = numpy.mean(self.values)
        Sig_y = numpy.std(self.values)

        # guess for lognormal
        s, scale = lognorm_to_scipyinput(M_y, Sig_y)
        m_x = math.log(scale)
        sig_x = s

        uni_b = 2*Sig_y
        uni_a = M_y - Sig_y

        return [[M_y, Sig_y],[m_x, sig_x], [uni_a,uni_b]]

    def make_cdfs(self):

        width = 1/len(self.values)
        y = []
        for i in range (len(self.values)):
            y.append(0+(width*i))

        cdf_list = [[stats.norm.cdf], [my_log_normal_cdf], [stats.uniform.cdf]]

        for i in range (len(cdf_list)):

            param = optimize.curve_fit(cdf_list[i][0], self.values, y, p0=self.init_guess[i])[0]

            cdf_list[i].append(param)


        cdf_list[1][1] = [cdf_list[1][1][1], 0, math.exp(cdf_list[1][1][0])]


        return cdf_list, y


    def ks_cdf(self):

        ks_list = []

        for i in range (len(self.dist_types)):
            ks = stats.kstest(self.values, self.dist_types[i], args=self.cdfs[i][1])[0]
            ks_list.append(ks)

        print('ks: ',ks_list)

        min = 2
        index = 0
        for i in range (len(ks_list)):
            if ks_list[i] < min:
                min = ks_list[i]
                index = i

        return index



#TODO get new version of scipy so that rv_histogram works, and then we can plot pdf, then fix make_pdf_hist

    def make_pdf_hist(self):

        #print(self.values)
        hist = numpy.histogram(self.values, bins=self.num_bins)
        hist_dist = stats.rv_histogram(hist)



        #bins, step = linspace(self.values[0],self.values[len(self.values)-1], self.num_bins, retstep=True)

        #plt.plot(bins,hist_dist.pdf(bins))
       #kde = stats.gaussian_kde(self.values)


        #return [bins, kde]


    def show_all_cdf(self):

        fig, ax = plt.subplots(1, 3, figsize=(12, 6))
        titles = ['CDF: With Normal Fit', 'CDF: With Log-Normal Fit', 'CDF: With Uniform Fit']
        totalwidth = 2*(max(self.values) - min(self.values))
        x1 = make_x1(totalwidth)
        print(x1)

        lognorm_param = self.cdfs[1][1]
        normal_param = self.cdfs[0][1]
        uniform_param = self.cdfs[2][1]
        print(lognorm_param,normal_param,uniform_param)



        for i in range (len(ax)):
            ax[i].set_xlabel('(ng/g) of ' + self.chem + ' in ' + self.animal, size='large')
            ax[i].set_ylabel('P(X < x)')
            ax[i].title.set_text(titles[i])
            ax[self.index].title.set_text(titles[self.index] + '(Considered Optimal by ks-test)')

        for i in range(len(ax)):
            ax[i].plot(self.values, self.y, 'ro', color='r')


        ax[0].plot(x1, stats.norm.cdf(x1, scale=normal_param[1], loc=normal_param[0]), linewidth=2.0, color='g')
        ax[1].plot(x1, stats.lognorm.cdf(x1, s=lognorm_param[0] ,scale=lognorm_param[2]), linewidth=2.0, color='g')
        ax[2].plot(x1, stats.uniform.cdf(x1, scale=uniform_param[1], loc=uniform_param[0]), linewidth=2.0, color='g')


        plt.show()


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
        #x = self.hist[0]

        if self.dist_types[self.index] == 'norm':
            mean = fit_params[0]
            std = fit_params[1]
            totalwidth = 2*(mean+(3*std))
            x1 = make_x1(totalwidth)

            #ax[0].set_xlim(xmin, xmax)
            #ax[0].plot(x, self.hist[1](x), 'ro', color='r')
            ax[0].plot(x1, stats.norm.pdf(x1, loc=mean, scale=std), linewidth=2.0)
            ax[1].set_xlim(self.values[0], self.values[len(self.values)-1])
            ax[1].plot(self.values,self.y, 'ro', color='r')
            ax[1].plot(x1, stats.norm.cdf(x1, loc=mean, scale=std), color='g')
            plt.show()

        if self.dist_types[self.index] == 'uniform':

            a = fit_params[0]
            b = fit_params[1]
            totalwidth = 2*(a+b)
            x1 = make_x1(totalwidth)

            #ax[0].set_xlim(xmin, xmax)
            #ax[0].plot(x, self.hist[1](x), 'ro', color='r')
            ax[0].plot(x1, stats.uniform.pdf(x1, loc=a, scale=b), linewidth=2.0)
            ax[1].set_xlim(self.values[0], self.values[len(self.values) - 1])
            ax[1].plot(self.values, self.y, 'ro', color='r')
            ax[1].plot(x1, stats.uniform.cdf(x1, loc=a, scale=b), color='g')
            plt.show()

        if self.dist_types[self.index] == 'lognorm':
            mu = fit_params[0]
            sig2 = fit_params[1]
            totalwidth = 2*(mu+(3*sig2))
            x1 = make_x1(totalwidth)

            #ax[0].set_xlim(xmin, xmax)
            #ax[0].plot(x, self.hist[1](x), 'ro', color='r')
            ax[0].plot(x1, stats.lognorm.pdf(x1, 1 ,loc=mu, scale=sig2), linewidth=2.0)
            ax[1].set_xlim(self.values[0], self.values[len(self.values) - 1])
            ax[1].plot(self.values, self.y, 'ro', color='r')
            ax[1].plot(x1, stats.lognorm.cdf(x1, 1, loc=mu, scale=sig2), color='g')
            plt.show()

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


###### graphing and lognormal help methods ##############


def make_x1(total_width):

    interval = total_width / float(1000)

    x1 = []
    for i in range(1000):
        toappend = interval * i
        x1.append(toappend)

    return x1

def my_log_normal_cdf(x, M_y, Sig_y):

    to_run = (numpy.log(x) - M_y)/Sig_y
    F_x = st.norm.cdf(to_run)

    return F_x


def lognorm_to_scipyinput(M_y,Sig_y):
    print(M_y,Sig_y)
    m_x = (2 * math.log(M_y)) - (.5) * (math.log(math.pow(Sig_y, 2) + math.pow(M_y, 2)))

    scale = math.exp(m_x)

    sigma2 = -2 * math.log(M_y) + math.log(math.pow(Sig_y, 2) + math.pow(M_y, 2))
    print(sigma2)
    s = math.sqrt(sigma2)

    return s, scale


#### overarching loop methods ##########

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

