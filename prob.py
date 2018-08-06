from scipy import stats as st
import pyDOE
import math

# TODO figure out how to avoid negative numbers when sampling from some distrubtions
# TODO Figure out the correct binning algorthim / fit cdf?

class loguniform:

    def __init__(self, loc=-1, scale=0, base=10):
        self.loc = loc
        self.scale = scale
        self.base = base

    def ppf(self,q):

        uniform = st.uniform(loc=self.loc, scale=self.scale)
        return math.pow(self.base, uniform.ppf(q))

class ResultDist:


    def __init__(self, values):

        self.values = values
        print(values)
        self.values.sort()
        self.bin_size = self.bin_size()
        self.hist = self.make_pdf_hist()

    def __str__(self):

        string = 'values: '
        for num in self.values:
            string += (str(num) + ' ')
        string += '\n\t\t\t'
        string += ('Bin Width: ' + str(self.bin_size))
        return string


    def bin_size(self):

        return (2*self.IQR())/(math.pow(len(self.values),(1/3)))

    def IQR(self):

        l = self.values[:len(self.values)//2]
        q1 = l[len(l)//2]
        l = self.values[len(self.values)//2:]
        q3 = l[len(l)//2]

        return q3-q1

    def make_pdf_hist(self):

        hist = {}
        num_bin = int((self.values[len(self.values) - 1] - self.values[0]) // self.bin_size) + 1

        for i in range(num_bin):
            hist[i] = 0

        print(self.bin_size)

        for i in range(len(self.values)):
            bottom = self.values[0]
            top = bottom + self.bin_size
            bin_key = 0
            print(self.values[i])
            for j in range(num_bin):
                if self.values[i] < top and self.values[i] >= bottom:
                    print('got in : ', end='')
                    hist[bin_key] += 1
                bottom = top
                top = bottom + self.bin_size
                bin_key += 1

        return hist




def set_hyper_cube(model_para, Var):

    v_iter = int(model_para[0])
    u_iter = int(model_para[1])
    bin_num = int(model_para[2])

    if Var.type == 'U':
        hype_sample = v_iter//bin_num
    else:
        hype_sample = u_iter//bin_num

    #print(hype_sample)

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
            mu = Var.param[0]
            math.exp(mu)
            sigma = Var.param[1]
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


def make_result_dist(dicts):

    result_dict = {}

    for region, values in dicts[0].items():
        result_dict[region] = {}
        for animal, values1 in values.items():
            result_dict[region][animal] = {}
            for chemical in values1.keys():
                values = get_values(region,chemical, animal, dicts)
                new_dist = ResultDist(values)
                result_dict[region][animal][chemical] = new_dist


def get_values(region,chem, animal, dicts):

    dist_list = []

    for dict in dicts:
        dist_list.append(dict[region][animal][chem])

    return dist_list

