from scipy import stats as st
from pyDOE import *
import Classes as cs

def set_hyper_cube(model_para, Var):

    v_iter = int(model_para[0])
    u_iter = int(model_para[1])
    bin_num = int(model_para[2])

    hype_sample_v = v_iter//bin_num
    hype_sample_u = u_iter//bin_num

    v_lhs = lhs(bin_num, samples=hype_sample_v)
    u_lhs = lhs(bin_num, samples=hype_sample_u)

    v_lhs = v_lhs.ravel()
    u_lhs = u_lhs.ravel()

    Var.u_lhs = u_lhs
    Var.v_lhs = v_lhs

def sample_dist(name, Var, i):

    if Var.dist == 'Normal':
        mean = Var.param[0]
        std = Var.param[1]
        if Var.type == 'U':
            point = Var.u_lhs[i]
        if Var.type == 'V':
            point = Var.v_lhs[i]
        sample = st.norm(loc=mean, scale=std).ppf(point)
        #sample = sampler.rvs(size=1)
    elif Var.dist == 'Uniform':
        a = Var.param[0]
        b = Var.param[1]
        if Var.type == 'U':
            point = Var.u_lhs[i]
        if Var.type == 'V':
            point = Var.v_lhs[i]
        sample = st.uniform(loc=a, scale=b).ppf(point)
    elif Var.dist == 'Triangle':
        a = Var.param[0]
        b = Var.param[1]
        c = (Var.param[2]-a)/b
        if Var.type == 'U':
            point = Var.u_lhs[i]
        if Var.type == 'V':
            point = Var.v_lhs[i]
        sample = st.triang(c, loc=a, scale=b).ppf(point)
    else:
        print(Var.dist)
        print('There is a unknown distribution in', '\''+name + '\'')
        exit(0)

    Var.value = float(sample)

