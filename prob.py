from scipy import stats as st
import Classes as cs

def sample_dist(Var, sample_num):

    if Var.dist == 'Normal':
        mean = Var.param[0]
        std = Var.param[1]
        samples = st.norm(loc=mean, scale=std).rvs(size=sample_num)
        return samples
    if Var.dist == 'Uniform':
        a = Var.param[0]
        b = Var.param[1]
        samples = st.uniform(loc=a, scale=b).rvs(size=sample_num)
        return samples
    if Var.dist == 'Triangle':
        a = Var.param[0]
        b = Var.param[1]
        c = (Var.param[2]-a)/b
        samples = st.triang(c, loc=a, scale=b).rvs(size=sample_num)
        return samples


def main():

    to_sample = cs.Var('U', 'Triangle', [0,1,.5])
    samples = sample_dist(to_sample, 400)
    #print(samples)

main()