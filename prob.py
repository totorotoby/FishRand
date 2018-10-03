import pyDOE
import math
import numpy
import matplotlib
from scipy import stats as st
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from scipy import optimize
from scipy import stats
from scipy import special

class Var:

    def __init__(self, type, dist, param):

        self.type = type   
        self.dist = dist   # distribution name
        self.param = param
        self.values = None
        self.lhs = None

    def __str__(self):

        to_print = '\ntype | ' + self.type + '\ndistribution | ' + self.dist + '\nparameters | ' + str(self.param)
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
            s, scale = lognorm_to_scipyinput(m_y, sig_y)
            self.values = st.lognorm(s=s, scale=scale).ppf(self.lhs)
        elif self.dist == 'Log-Uniform':
            a = self.param[0]
            b = self.param[1]
            self.values = Loguniform(loc=a, scale=b).ppf(self.lhs)
        elif self.dist == 'Beta':
            alpha = self.param[0]
            beta = self.param[1]
            self.values = st.beta(alpha, beta).ppf(self.lhs)
            self.values = sorted(self.values)
        elif self.dist == 'Weibull':
            lamb = self.param[0]
            k = self.param[1]
            self.values = st.weibull_min(c=k, scale=lamb).ppf(self.lhs)
        else:
            print(self.dist)
            print('There is a unknown distribution called ', '\'' + self.dist + '\'')
            exit(0)


class Loguniform:

    def __init__(self, loc=-1, scale=0, base=math.e):
        self.loc = loc
        self.scale = scale
        self.base = base

    def ppf(self, q):

        uniform = st.uniform(loc=self.loc, scale=self.scale)
        return math.pow(self.base, uniform.ppf(q))


class ResultDist:

    dist_types = ['norm', 'lognorm', 'uniform', 'gamma']

    def __init__(self, values, chemical, animal):
        self.count = 0

        self.chem = chemical
        self.animal = animal
        self.values = values
        self.v_mean_std = [round(numpy.mean(values), 4), round(numpy.std(values), 4)]
        self.display = 4
        self.values.sort()
        self.init_guess = self.inital_guess()
        self.num_bins = len(self.values)//50
        self.hist = self.make_pdf_hist()
        self.cdfs, self.y = self.make_cdfs()
        self.pdfs = [stats.norm.pdf, stats.lognorm.pdf, stats.uniform.pdf, stats.gamma.pdf]
        self.cdf_list = [stats.norm.cdf, stats.lognorm.cdf, stats.uniform.cdf, stats.gamma.cdf]
        self.x1 = make_x1(self.values)
        self.index = self.ks_cdf()
        self.best_para = self.bestparam()

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

        gamma_k = (math.pow(M_y, 2))/(math.pow(Sig_y, 2))
        gamma_theta = (math.pow(Sig_y, 2))/M_y
        
        # returns guess array...order: [normal, lognormal, uniform, gamma]
        return [[M_y, Sig_y], [m_x, sig_x], [uni_a, uni_b], [gamma_k, gamma_theta]]

    def make_cdfs(self):

        width = 1/len(self.values)
        y = []
        for i in range(len(self.values)):
            y.append(0+(width*i))

        cdf_list = [[stats.norm.cdf], [my_log_normal_cdf], [stats.uniform.cdf], [my_gamma_cdf]]

        for i in range(len(cdf_list)):

            param = optimize.curve_fit(cdf_list[i][0], self.values, y, p0=self.init_guess[i])[0]
            cdf_list[i].append(param)

        cdf_list[1][1] = [cdf_list[1][1][1], 0, math.exp(cdf_list[1][1][0])]
        cdf_list[3][1] = [cdf_list[3][1][0], 0, cdf_list[3][1][1]]

        return cdf_list, y

    def ks_cdf(self):

        ks_list = []

        for i in range(len(self.dist_types)):
            ks = stats.kstest(self.values, self.dist_types[i], args=self.cdfs[i][1])[0]
            ks_list.append(ks)

        minimum = 2
        index = 0
        for i in range(len(ks_list)):
            if ks_list[i] < minimum:
                minimum = ks_list[i]
                index = i

        return index

    def make_pdf_hist(self):

        hist, bins = numpy.histogram(self.values, bins=self.num_bins, normed=True)
        return [hist, bins]
     
    def plot_cdf(self):

        fig, ax = plt.subplots(1, 4, figsize=(24, 12))

        titles = ['CDF: With Normal Fit', 'CDF: With Log-Normal Fit', 'CDF: With Uniform Fit',
                  'CDF: With Gamma Fit']

        # CDF x-axis
        x1 = self.x1
       
        # CDF parameters
        lognorm_param = self.cdfs[1][1]
        normal_param = self.cdfs[0][1]
        uniform_param = self.cdfs[2][1]
        gamma_param = self.cdfs[3][1]

        for i in range(len(ax)):
            
            # CDF axis labels
            ax[i].set_xlabel('(ng/g) of ' + self.chem + ' in ' + self.animal, size='large')
            ax[i].set_ylabel('P(X < x)')

            # titles
            ax[i].title.set_text(titles[i])
            ax[self.index].title.set_text(titles[self.index] + '(Considered Optimal by ks-test)')

            # plot points
            ax[i].plot(self.values, self.y, 'ro', color='r')

        ax[0].plot(x1, stats.norm.cdf(x1, scale=normal_param[1], loc=normal_param[0]), linewidth=2.0, color='g')
        ax[1].plot(x1, stats.lognorm.cdf(x1, s=lognorm_param[0], scale=lognorm_param[2]), linewidth=2.0, color='g')
        ax[2].plot(x1, stats.uniform.cdf(x1, scale=uniform_param[1], loc=uniform_param[0]), linewidth=2.0, color='g')
        ax[3].plot(x1, stats.gamma.cdf(x1, a=gamma_param[0], loc=gamma_param[1], scale=gamma_param[2]),
                   linewidth=2.0, color='g')

        fig.tight_layout()
        
    def plot_pdf(self):

        titles = ['PDF: With Normal Fit', 'PDF: With Log-Normal Fit', 'PDF: With Uniform Fit',
                  'PDF: With Gamma Fit']

        fig1, ax1 = plt.subplots(1, 4, figsize=(24, 12))

        x1 = self.x1

        lognorm_param = self.cdfs[1][1]
        normal_param = self.cdfs[0][1]
        uniform_param = self.cdfs[2][1]
        gamma_param = self.cdfs[3][1]

        for i in range(len(ax1)):
            ax1[i].set_xlabel('(ng/g) of ' + self.chem + ' in ' + self.animal, size='large')
            ax1[i].set_ylabel('P(x)')
            ax1[self.index].title.set_text(titles[self.index] + '(Considered Optimal by ks-test)')
            ax1[i].title.set_text(titles[i])
            ax1[i].scatter(self.hist[1][:-1], self.hist[0], s=16)
            ax1[i].set_ylim(min(self.hist[0]))

        ax1[0].plot(x1, stats.norm.pdf(x1, scale=normal_param[1], loc=normal_param[0]), linewidth=2.0, color='g')
        ax1[1].plot(x1, stats.lognorm.pdf(x1, s=lognorm_param[0], scale=lognorm_param[2]), linewidth=2.0, color='g')
        ax1[2].plot(x1, stats.uniform.pdf(x1, scale=uniform_param[1], loc=uniform_param[0]), linewidth=2.0, color='g')
        ax1[3].plot(x1, stats.gamma.pdf(x1, a=gamma_param[0], loc=gamma_param[1], scale=gamma_param[2]),
                    linewidth=2.0, color='g')

        fig1.tight_layout()

    def plot_single(self, temp_index):

        titles = ['CDF: With Normal Fit', 'CDF: With Log-Normal Fit', 'CDF: With Uniform Fit', 'CDF: With Gamma Fit']
        titles1 = ['PDF: With Normal Fit', 'PDF: With Log-Normal Fit', 'PDF: With Uniform Fit', 'PDF: With Gamma Fit']

        fig, ax = plt.subplots(1, 2, figsize=(24, 12))

        x1 = self.x1
        
        temp_param = self.cdfs[temp_index][1]

        for i in range(len(ax)):
            ax[i].set_xlabel('(ng/g) of ' + self.chem + ' in ' + self.animal, size='large')
            
        ax[0].title.set_text(titles1[temp_index])
        ax[0].set_ylabel('P(x)')
        ax[0].scatter(self.hist[1][:-1], self.hist[0], s=16)
        ax[0].set_ylim(min(self.hist[0]))
        if len(temp_param) == 2:
            ax[0].plot(x1, self.pdfs[temp_index](x1, scale=temp_param[1], loc=temp_param[0]), linewidth=2.0, color='g')
        if len(temp_param) == 3:
            ax[0].plot(x1, self.pdfs[temp_index](x1, temp_param[0], loc=temp_param[1], scale=temp_param[2]),
                       linewidth=2.0, color='g')

        ax[1].title.set_text(titles[temp_index])
        ax[1].set_ylabel('P(X < x)')
        ax[1].plot(self.values, self.y, 'ro', color='r')
        if len(temp_param) == 2:
            ax[1].plot(x1, self.cdf_list[temp_index](x1, scale=temp_param[1], loc=temp_param[0]), linewidth=2.0, color='g')
        if len(temp_param) == 3:
            ax[1].plot(x1, self.cdf_list[temp_index](x1, temp_param[0], loc=temp_param[1], scale=temp_param[2]),
                       linewidth=2.0, color='g')

        fig.tight_layout()

    def show(self, temp_index):

        if self.display == 0:
            self.plot_cdf()
        elif self.display == 1:
            self.plot_pdf()
        elif self.display == 2:
            self.plot_pdf()
            self.plot_cdf()
        else:
            self.plot_single(temp_index)
        plt.show()

    def bestparam(self):

        params = []
        self.count += 1
        string = ''
        string += str(self.dist_types[self.index] + '(')

        for i in range(len(self.cdfs[self.index][1])-1):
            params.append(self.cdfs[self.index][1][i])
            string += str("{0:.4f}".format(self.cdfs[self.index][1][i]))
            string += ', '
        params.append(self.cdfs[self.index][1][i+1])
        string += str("{0:.4f}".format(self.cdfs[self.index][1][i+1]))
        string += ')'
        return [string, params]


# graphing, gamma, and lognormal help methods
def make_x1(values):

    mean = numpy.mean(values)
    std = numpy.std(values)
    
    x1 = numpy.linspace(mean - (4*std), mean + (6*std), num=200)

    return x1

def my_log_normal_cdf(x, M_y, Sig_y):

    to_run = (numpy.log(x) - M_y)/Sig_y
    F_x = st.norm.cdf(to_run)

    return F_x


def my_gamma_cdf(x, alpha, beta):

    bottom = special.gamma(alpha)
    top = special.gammainc(alpha, beta*x)

    return top/bottom


def lognorm_to_scipyinput(M_y,Sig_y):

    m_x = (2 * math.log(M_y)) - (.5) * (math.log(math.pow(Sig_y, 2) + math.pow(M_y, 2)))

    scale = math.exp(m_x)

    sigma2 = -2 * math.log(M_y) + math.log(math.pow(Sig_y, 2) + math.pow(M_y, 2))
    try:
        s = math.sqrt(sigma2)
    except:
        print('normal to log normal bug', sigma2, M_y, Sig_y)
        exit(0)

    return s, scale


# overarching loop methods

def set_hyper_samp_cube(model_para, Var):

    u_iter = int(model_para[0])
    v_iter = int(model_para[1])
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
                values = get_values(region, chemical, animal, dicts)
                new_dist = ResultDist(values, chemical, animal)
                result_dict[region][animal][chemical] = new_dist

    return result_dict


def get_values(region, chem, animal, dicts):

    dist_list = []

    for diction in dicts:
        dist_list.append(diction[region][animal][chem])

    return dist_list
