# Where we run the bioaccumlation model
import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import Classes as obj
import FR_Input_Output
import prob as pr

def set_all_h_and_s(model_para, all_data):

    for list in all_data:
        for item in list:
            for i in range(len(item)):
                if type(item[i]) == pr.Var:
                    pr.set_hyper_samp_cube(model_para, item[i])

def check_inst_non_st(inst,u_iter, v_iter):

    new_region = []
    for j in range (len(inst)):
        if type(inst[j]) == pr.Var:
            if inst[j].type == 'U':
                new_region.append(inst[j].values[u_iter])
            if inst[j].type == 'V':
                new_region.append(inst[j].values[v_iter])
        else:
            new_region.append(inst[j])

    return new_region


# initiates all regions
def init_region(reg_data,u_iter, v_iter):

    regions = []

    for i in range(len(reg_data)):
        region = reg_data[i]
        region = check_inst_non_st(region,u_iter, v_iter)
        toadd = obj.Region(region[0],region[1],region[2],region[3],region[4],region[5],region[6])

        # dealing with cox
        if region[7] != '':
            toadd.set_cox(region[7])
        else:
            toadd.calc_cox()

        # dealing with adoc apoc
        if region[8] and region[9] != '':
            toadd.set_adoc_apoc(region[8],region[9])

        # checking to see if anything is missing
        if toadd.init_check():
            regions.append(toadd)
        else:
            print("Something is wrong with your ", i, " Regional Entry")
            exit(0)

    return regions

# initiates chemicals for a specific single region
def init_chems(chem_data, region,u_iter, v_iter):

    chemicals = []
    for i in range(len(chem_data)):

        chemical = chem_data[i]
        chemical = check_inst_non_st(chemical,u_iter, v_iter)

        toadd = obj.Chemical(chemical[0],chemical[1],chemical[2])

        # dealing with cwto and cwdo
        if chemical[3] != '':
            toadd.set_cwto(chemical[3])
        if chemical[4] != '':
            toadd.set_cwdo(chemical[4])

        # dealing with ddoc and dpoc
        if chemical[5] and chemical[6] != '':
            toadd.set_ddoc_dpoc(chemical[5],chemical[6])

        # calculating Phi
        toadd.calc_phi_and_cwdo(region)
        toadd.calc_pore_water(region)

        # checking to see if anything is missing
        if toadd.init_check():
            chemicals.append(toadd)
        else:
            print("Something is wrong with your ",i ," Chemical Entry")
            exit(0)


    return chemicals


def init_phyto(phyto_data, chemicals,u_iter, v_iter):

    phytos = []
    phyto = phyto_data[0]
    phyto = check_inst_non_st(phyto,u_iter, v_iter)

    # Potential loop in the future

    toadd = obj.Pplank(phyto[0],phyto[1])

    if phyto[2] != '':
        toadd.set_vlp(phyto[2])
    if phyto[3] != '':
        toadd.set_vnp(phyto[3])

    toadd.calc_vwp()

    for chemical in chemicals:

        kow = chemical.Kow
        k_1 = toadd.calc_k1(kow)
        toadd.k_1.append(k_1)
        k_2 = toadd.calc_k2(kow,k_1)
        toadd.k_2.append(k_2)

    if toadd.init_check():
        phytos.append(toadd)
    else:
        print("Something is wrong with your ", 1,  "st Phyto Entry")
        exit(0)


    return phytos


def init_zoop(zoo_data, region, phyto, chemicals,u_iter, v_iter):

    zoops = []
    zoop = zoo_data[0]
    zoop = check_inst_non_st(zoop,u_iter, v_iter)

    # Potential loop in the future

    toadd = obj.Zooplank(zoop[0],zoop[1],zoop[2], zoop[10])

    # inital set
    if zoop[3] != '':
        toadd.set_vnb(zoop[3])
    if zoop[4] != '':
        toadd.set_mp(zoop[4])
    if zoop[5] != '' and zoop[6] != '' and zoop[7] != '':
        toadd.set_el_en_ew(zoop[5],zoop[6],zoop[7])
    if zoop[8] != '':
        toadd.set_gd(zoop[8])
    if zoop[9] != '':
        toadd.set_kg(zoop[9])

    # second check
    toadd.calc_mo()
    toadd.calc_vwb()
    if toadd.Gd == None:
        toadd.calc_gd(region.T,region.Css, 0)
    if toadd.Kg == None:
        toadd.calc_kg(region.T)

    toadd.calc_diet_per(phyto)
    toadd.calc_gut_per()
    toadd.calc_gv(region)
    toadd.calc_gf()

    for chemical in chemicals:

        kow = chemical.Kow
        ed = chemical.Ed
        ew = chemical.Ew
        k_1 = toadd.calc_k1(ew)
        toadd.k_1.append(k_1)
        k_2 = toadd.calc_k2(kow, k_1)
        toadd.k_2.append(k_2)
        k_gb = toadd.calc_kgb(kow)
        toadd.k_gb.append(k_gb)
        k_e = toadd.calc_ke(ed,k_gb)
        toadd.k_e.append(k_e)
        k_d = toadd.calc_kd(ed)
        toadd.k_d.append(k_d)


    if toadd.init_check():
        zoops.append(toadd)
    else:
        print("Something is wrong with your ", 1,  "st Zoop Entry")
        exit(0)

    return zoops


def init_fish(fish_data, diet_data, region, chemicals, phyto, zoop,u_iter, v_iter):

    tempfishs = []
    tempfishs.append(zoop)
    tempfishs.append(phyto)
    fishs = []


    for fish in fish_data:

        fish = check_inst_non_st(fish,u_iter, v_iter)

        diet = diet_data[fish[0]]

        toadd = obj.Fish(fish[0],fish[1],fish[2],fish[4], diet, fish[10])
        tempadd = obj.Fish(fish[0],fish[1],fish[2],fish[4], diet, fish[10])
        # initial variables
        if fish[3] != '':
            toadd.set_vnb(fish[3])
            tempadd.set_vnb(fish[3])

        tempadd.calc_vwb()
        tempfishs.append(tempadd)

        if fish[5] != '':
            toadd.set_gd(fish[5])
        if fish[6] != '' and fish[7] != '' and fish[8] != '':
            toadd.set_el_en_ew(fish[6],fish[7],fish[8])
        if fish[9] != '':
            toadd.set_kg(fish[9])

        # second variables
        toadd.calc_mo()
        toadd.calc_vwb()
        toadd.calc_gv(region)
        if toadd.Gd == None:
            toadd.calc_gd(region.T, region.Css, toadd.flag)
        if toadd.Kg == None:
            toadd.calc_kg(region.T)
        fishs.append(toadd)

    for i in range(len(fishs)):
        fishs[i].calc_diet_per(tempfishs,region)
        fishs[i].calc_gut_per()
        fishs[i].calc_gf()

    count = 0
    # going back through so that we can now set diets
    for i in range (len(fishs)):
        fish = fishs[i]
        k_1 = []
        k_2 = []
        k_gb = []
        k_e = []
        k_d = []

        for j in range (len(chemicals)):
            kow = chemicals[j].Kow
            ew = chemicals[j].Ew
            ed = chemicals[j].Ed
            k1 =fish.calc_k1(ew)
            k_1.append(k1)
            k2 = fish.calc_k2(kow,k1)
            k_2.append(k2)
            kgb = fish.calc_kgb(kow)
            k_gb.append(kgb)
            ke = fish.calc_ke(ed,kgb)
            k_e.append(ke)
            kd = fish.calc_kd(ed)
            k_d.append(kd)


        fish.k_1 = k_1
        fish.k_2 = k_2
        fish.k_gb = k_gb
        fish.k_e = k_e
        fish.k_d = k_d


        if fish.init_check():
            count += 1

    if count == len(fishs):
        return fishs
    else:
        print('Something is wrong with a fish entry')
        exit(0)


def reorder_fish(fishs):

    new_order = []
    below = []
    # find the root fish
    while True:
        fish = find_base_fish(fishs)
        if type(fish) == obj.Fish:
            new_order.append(fish)
            below.append(fish.name)
            fishs.remove(fish)
        if fish == None:
            break

        below.append('Phytoplankton')
        below.append('Zooplankton')
        below.append('Sediment/Detritus')

    while True:
        nextfish = find_next_fish(fishs, new_order, below)
        if type(nextfish) == obj.Fish:
            new_order.append(nextfish)
            below.append(nextfish.name)
            fishs.remove(nextfish)
        if nextfish == None:
            break

    return new_order


def find_base_fish(fishs):

    basefish = None

    for fish in fishs:
        # counter for non 0 non phyto zoop diet
        count = 0
        for entry in fish.diet_frac:

            if (entry[0] != 'Sediment/Detritus' and entry[0] != 'Phytoplankton' and entry[0] != 'Zooplankton') and entry[1] > 0:
                count += 1

        if count == 0:
            basefish = fish

    return basefish


def find_next_fish(possible_fish, current_fish, below):

    nextfish = None

    for fish in possible_fish:
        count = 0
        for entry in fish.diet_frac:
            if (entry[0] not in below) and entry[1] > 0:

                count += 1

        if count == 0:
            nextfish = fish


    return nextfish


def solve(regions, chemicals, phytos, zoops, fishs):

    # nested dictonary where can look up first by region then chemical then animal to find concentration
    #assuming one region for now
    conc_log = {}
    conc_log[regions[0].name] = {}

    conc_log = solve_phyto(chemicals,phytos,conc_log, regions[0])
    conc_log = solve_zoop(chemicals,zoops, conc_log, regions[0])
    conc_log = solve_fish(chemicals,fishs,conc_log,regions[0])

    return conc_log



def solve_phyto(chemicals,phytos, conc_log, region):

    phytolog = conc_log[region.name]
    Cb = []
    #assuming one phyto
    phyto = phytos[0]
    phytolog[phyto.name] = {}
    for i in range (len(chemicals)):
        cwd = chemicals[i].Cwdo
        conc_in_phyto = phyto.solve_steady_state(cwd, i)
        Cb.append(conc_in_phyto)
        phytolog[phyto.name][chemicals[i].name] = conc_in_phyto
    phyto.Cb = Cb

    return conc_log


def solve_zoop(chemicals, zoops, conc_log, region):

    zooplog = conc_log[region.name]
    #assuming one zoop
    zoops = zoops[0]
    Cb = []
    zooplog[zoops.name] = {}
    for i in range (len(chemicals)):

        phi = chemicals[i].phi
        Cwto = chemicals[i].Cwto
        Cwds = chemicals[i].Cwdo
        phyto_con = conc_log[region.name]['Phytoplankton'][chemicals[i].name]

        conc_in_zoops = zoops.solve_steady_state(phi, i, Cwto, Cwds, phyto_con)
        Cb.append(conc_in_zoops)
        zooplog[zoops.name][chemicals[i].name] = conc_in_zoops

    zoops.Cb = Cb

    return conc_log



def solve_fish(chemicals, fishs, conc_log, region):

    fishlog = conc_log[region.name]

    for i in range (len(fishs)):
        fishlog[fishs[i].name] = {}
        Cb = []
        for j in range (len(chemicals)):
            phi = chemicals[j].phi
            Cwp = chemicals[j].Cwp
            Cwds = chemicals[j].Cwdo
            con_in_i = fishs[i].solve_steady_state(phi, j, Cwp, Cwds, fishlog, chemicals[j])
            Cb.append(con_in_i)
            fishlog[fishs[i].name][chemicals[j].name] = con_in_i
        fishs[i].Cb = Cb

    return conc_log


def pretty(d, indent=0):
   for key, value in d.items():
      print('\t' * indent + str(key))
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print('\t' * (indent+1) + str(value))



def run_bio(flag, filename, endname):

    if flag == 0:
        all_data = FR_Input_Output.convert_to_lists(filename)[1]
        conc_log = single_iter(all_data[0], all_data[1], all_data[2], all_data[3], all_data[4], all_data[5],0, endname)
        return conc_log
    else:
        dictionares = []
        model_para, all_data  = FR_Input_Output.convert_to_lists(filename)

        v_iter = int(model_para[0])
        u_iter = int(model_para[1])

        set_all_h_and_s(model_para, all_data)
        inner_count = 0
        u_count = 0
        print('percentage done: ')
        while (u_count < u_iter):
            u_count += 1
            v_count = 0
            while (v_count < v_iter):
                log = single_iter(all_data[0], all_data[1], all_data[2], all_data[3], all_data[4], all_data[5],1 ,endname ,u_iter=u_count,v_iter=inner_count)
                dictionares.append(log)
                v_count += 1
                inner_count += 1
                print( '\r' + str(100*(inner_count/(v_iter*u_iter))), end='')



        results_dic = pr.make_result_dist(dictionares)

        return results_dic

def result_print(results_dic):

    for region in results_dic.values():
        for animal, animals in region.items():
            for chemical, cont in animals.items():
                print(cont.bestparam()[0])
                cont.plot_info()


def plot_small(graph_list):

    width = 1 / len(graph_list)
    y = []
    for i in range(len(graph_list)):
        y.append(0 + (width * i))
    print(len(graph_list[0]),len(graph_list))
    for j in range(len(graph_list[0])):
        x = []
        for i in range(len(graph_list)):
            x.append(graph_list[i][j])
        x = sorted(x)
        print(x)
        plt.plot(x,y,'ro', color='r')
        plt.xlim(x[0],x[len(x)-1])
        plt.show()

def single_iter(reg_data, chem_data, fish_data, zoo_data, phyto_data, diet_data,flag, inputfilename, u_iter=0, v_iter=0):

    regions = init_region(reg_data,u_iter, v_iter)
    chemicals = init_chems(chem_data, regions[0],u_iter, v_iter)
    phytos = init_phyto(phyto_data, chemicals,u_iter, v_iter)
    zoops = init_zoop(zoo_data, regions[0], phytos[0], chemicals,u_iter, v_iter)
    fishs = init_fish(fish_data, diet_data, regions[0], chemicals, phytos[0], zoops[0],u_iter, v_iter)
    fishs = reorder_fish(fishs)
    conc_log = solve(regions, chemicals, phytos, zoops, fishs)
    if flag == 0:
        anwser = str(input('\nWould you like to save results to an excel sheet?\n'))
        if anwser == 'y':
            FR_Input_Output.deter_write_output(regions, fishs, chemicals, phytos, zoops, inputfilename)
        return conc_log
    else:
        return conc_log
