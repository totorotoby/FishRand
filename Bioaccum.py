# Where we run the bioaccumlation model

import Classes as obj
import FR_Input


# initiates all regions
def init_region(reg_data):

    regions = []

    for i in range(len(reg_data)):
        region = reg_data[i]
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
def init_chems(chem_data, region):

    chemicals = []
    for i in range(len(chem_data)):
        chemical = chem_data[i]
        toadd = obj.Chemical(chemical[0],chemical[1],chemical[2])

        # dealing with cwto and cwdo
        if chemical[3] and chemical[4] != '':
            toadd.set_cwto_cwdo(chemical[3],chemical[4])

        # dealing with ddoc and dpoc
        if chemical[5] and chemical[6] != '':
            toadd.set_ddoc_dpoc(chemical[5],chemical[6])

        # calculating Phi
        toadd.calc_phi(region)

        # checking to see if anything is missing
        if toadd.init_check():
            chemicals.append(toadd)
        else:
            print("Something is wrong with your ",i ," Chemical Entry")
            exit(0)


    return chemicals


def init_phyto(phyto_data, chemicals):

    phytos = []
    phyto = phyto_data

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


def init_zoop(zoo_data, region, phyto, chemicals):

    zoops = []
    zoop = zoo_data

    # Potential loop in the future

    toadd = obj.Zooplank(zoop[0],zoop[1],zoop[2])

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
        k_1 = toadd.calc_k1(kow)
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


def init_fish(fish_data, diet_data, region, chemicals, phyto, zoop):

    tempfishs = []
    tempfishs.append(zoop)
    fishs = []


    for fish in fish_data:

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

        if toadd.Gd == None:
            toadd.calc_gd(region.T, region.Css)
        if toadd.Kg == None:
            toadd.calc_kg(region.T)
        fishs.append(toadd)

    for i in range(len(fishs)):
        fishs[i].calc_diet_per(tempfishs)
        fishs[i].calc_gut_per()
        fishs[i].calc_gv(region)
        fishs[i].calc_gf()

    # going back through so that we can now set diets
    count = 0
    # bugged ! #
    for i in range(len(fishs)):

        for chemical in chemicals:

            kow = chemical.Kow
            ed = chemical.Ed
            k_1 = fishs[i].calc_k1(kow)
            fishs[i].k_1.append(k_1)
            k_2 = fishs[i].calc_k2(kow, k_1)
            fishs[i].k_2.append(k_2)
            k_gb = fishs[i].calc_kgb(kow)
            fishs[i].k_gb.append(k_gb)
            k_e = fishs[i].calc_ke(ed, k_gb)
            fishs[i].k_e.append(k_e)
            k_d = fishs[i].calc_kd(ed)
            fishs[i].k_d.append(k_d)

    # bugged ! #


        if fishs[i].init_check():
            count += 1

    print(len(fishs[0].k_1))
    if count == len(fishs):
        return fishs
    else:
        print('Something is wrong with a fish entry')
        exit(0)




def solve(regions, chemcicals, phytos, zoops, fishs):

    # nested dictonary where can look up first by chemical then animal to find concentration
    conc_log = {{}}



def main():
    reg_data, chem_data, fish_data, zoo_data, phyto_data, diet_data = FR_Input.convert_to_lists("FR_Input.xlsx")

    regions = init_region(reg_data)

    chemicals = init_chems(chem_data,regions[0])

    phytos = init_phyto(phyto_data,chemicals)

    zoops = init_zoop(zoo_data,regions[0],phytos[0],chemicals)

    fishs = init_fish(fish_data, diet_data, regions[0],chemicals,phytos[0],zoops[0])


main()