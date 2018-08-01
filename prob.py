from scipy import stats as st
import pyDOE
import Classes as cs

# TODO figure out how to avoid negative numbers when sampling from some distrubtions
# TODO Figure out the correct binning algorthim / fit cdf?


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
        # TODO add more distributions
        if Var.dist == 'Normal':
            mean = Var.param[0]
            std = Var.param[1]
            point = Var.lhs[i]
            sample = st.norm(loc=mean, scale=std).ppf(point)
            #sample = sampler.rvs(size=1)
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
        else:
            print(Var.dist)
            print('There is a unknown distribution in', '\''+name + '\'')
            exit(0)


        Var.value = float(sample)

