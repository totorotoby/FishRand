# The fish, invert, Pplank, Zooplank, and chemical objects has parameters for input that are from Arnot (2004) right now
import math
#
class Fish:

    def __init__(self, name, weight, vlb, mp, diet_frac, flag, times, num_regions, num_chemicals, vnb=.2, e_l=.75, e_n=.75, e_w=.5):
        self.name = name
        self.Wb = weight  # wet weight
        self.Vlb = vlb  # Percent Lipid Content
        self.Vnb = vnb  # percent Nonlipid organic matter
        self.Mp = mp  # Percent Water Ventilated
        self.diet_frac = diet_frac
        self.e_l = e_l  # Dietary absorption efficiency of lipids
        self.e_n = e_n  # Dietary absorption efficiency of nonlipid organic matter
        self.e_w = e_w  # Dietary absorption efficiency of water
        self.flag = flag
        self.Mo = 0  # Percent Over Water Ventilated
        self.Vwb = 0  # Percent Water Content
        self.Vld = 0  # lipid fraction of diet
        self.Vndc = [0 for _ in range(num_regions)]  # nonlip fraction of diet from dirt
        self.Vndm = 0   # nonlip fracion of diet from animal
        self.Vwd = 0  # Water fraction of diet
        self.Vlg = 0  # lipid fraction of gut
        self.Vng = 0  # Non-lipid fraction of gut
        self.Vwg = 0  # water fraction of gut

        self.gd_set = 0
        self.kg_set = 0

        self.Gv = [[0 for _ in range(num_regions)] for _ in range(times)]  # Gill ventilation rate
        self.Gd = [[0 for _ in range(num_regions)] for _ in range(times)]  # Feeding rate
        self.Kg = [[0 for _ in range(num_regions)] for _ in range(times)]  # Growth rate
        self.Gf = [[0 for _ in range(num_regions)] for _ in range(times)]  # fecal egestion rate
        self.k_1 = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in
                    range(times)]  # clearance rate constant
        self.k_2 = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in
                    range(times)]  # (rate constant chemical elem via respiratory)
        self.k_gb = [0 for _ in range(num_chemicals)]  # Gut biota partition coefficient
        self.k_e = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in
                    range(times)]  # (rate constant via excretion into egested feces)
        self.k_d = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in
                    range(times)]  # (clearance rate constant via ingestion of food)

        self.Mb = 0  # mass of chemical for multiple times
        self.Cb = 0  # concentration for steady state

    def __str__(self):
        return self.__dict__.__str__()


    def set_el_en_ew(self, el, en, ew):
        self.e_l = el
        self.e_n = en
        self.e_w = ew

    def set_gd(self,gd, times, len_reg):
        for i in range (times):
            for j in range(len_reg):
                self.Gd[i][j] = gd
        self.gd_set = 1

    def set_kg(self, kg, times, len_reg):

        for i in range (times):
            for j in range(len_reg):
                self.Kg[i][j] = kg
        self.kg_set = 1

    def set_vnb(self,vnb):
        self.Vnb = vnb

    def calc_vwb(self):
        self.Vwb = (1 - (self.Vlb + self.Vnb))

    def calc_mo(self):
        self.Mo = (1 - self.Mp)

    # sets the Gill ventilation rate of zooplank
    def calc_gv(self, regional_cox, time, region_index):
        try:
            self.Gv[time][region_index] = 1400 * ( math.pow(self.Wb, .65) / regional_cox[time])
        except ValueError:
            print('Weight of ' + self.name + ' was sampled as a negative number. This makes no sense. Check the range of the weight distribution. \n Caught at Gv calculation.')

    # sets G_d for this zooplank
    def calc_gd_filter(self, css, time_index ,region_index, sigma=1):
        self.Gd[time_index][region_index] = self.Gv[time_index][region_index] * css * sigma

    def calc_gd_no_filter(self, T, time_index, region_index):
        temp_at_time = T[time_index]
        self.Gd[time_index][region_index] = .022 * math.pow(self.Wb, .85) * math.exp((.06 * temp_at_time))

    def calc_kg(self, T_region, time_index, region_index):
        if T_region[time_index] < 15:
            self.Kg[time_index][region_index] = .0005 * math.pow(self.Wb, -.2)
        if T_region[time_index] >= 15:
            self.Kg[time_index][region_index] = .00251 * math.pow(self.Wb, -.2)

    # sets Gf for this zooplank
    def calc_gf(self):
        self.Gf = (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc) + ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd)) * self.Gd

    # returns k_1 of chemical for this zooplank
    def calc_k1(self, chem_ew, time_index, region_index, chemical_index):
        self.k_1[time_index][region_index][chemical_index] = chem_ew * (self.Gv[time_index][region_index] / self.Wb)

    # returns k_2 of chemical for this zooplank
    def calc_k2(self, chem_kow, k_1, time_index, region_index, chemical_index, beta=.035, density_lip=.9 , density_w=1):
        k_bw = (self.Vlb * chem_kow)/density_lip + (self.Vnb * beta * chem_kow) + (self.Vwb/density_w)
        self.k_2[time_index][region_index][chemical_index] = k_1[time_index][region_index][chemical_index] / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, ed, time_index, region_index, chemical_index):

        self.k_d[time_index][region_index][chemical_index] = ed * (self.Gd[time_index][region_index] / self.Wb)

    # sets the percentages of the gut
    def calc_gut_per(self):

        self.Vlg = ((1 - self.e_l) * self.Vld) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vngc = ((1 - self.e_n) * self.Vndc) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vngm = ((1 - self.e_n) * self.Vndm) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vwg = ((1 - self.e_w) * self.Vwd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))





    # returns K_gb for certian chemical
    def calc_kgb(self, chem_kow, chem_index,beta1=.35, beta2=.035, density_lip=.9, density_w=1, z_water=.05):
        top = ((self.Vlg * (z_water*chem_kow))/density_lip) + (self.Vng * beta1 * (z_water*chem_kow)) + (z_water*self.Vwg/density_w)
        bottom = ((self.Vlb * (z_water*chem_kow))/density_lip) + (self.Vnb * beta2 * z_water * chem_kow) + (z_water*self.Vwb/density_w)
        self.k_gb[chem_index] = top/bottom


    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, time_index, region_index,  chem_index):

        self.k_e[time_index][region_index][chem_index] = (self.k_gb[chem_index]/self.Wb) * chem_ed * self.Gf[time_index][region_index]


    # sets the percentages of zooplank diet that are lipid, non-lipid and water
    def calc_diet_per(self, fishlog, Ocs, region_index, region_len):


        total_nonlip_c = 0
        total_nonlip_m = 0
        total_lip = 0
        count = 0

        # loop through possible types of fishes
        for i in range(len(fishlog)):
            # loop through fish names in diet
            for j in range(len(self.diet_frac)):
                # finding right fish

                if fishlog[i].name == self.diet_frac[j][0] and self.diet_frac[j][0] == 'Phytoplankton':

                    l_toadd = self.diet_frac[j][1] * (fishlog[i].Vlb)
                    nl_toadd = self.diet_frac[j][1] * (fishlog[i].Vnb)

                    total_nonlip_c += nl_toadd
                    total_lip += l_toadd


                elif fishlog[i].name == self.diet_frac[j][0]:

                    l_toadd = self.diet_frac[j][1] * (fishlog[i].Vlb)
                    nl_toadd = self.diet_frac[j][1] * (fishlog[i].Vnb)

                    total_nonlip_m += nl_toadd
                    total_lip += l_toadd

                elif self.diet_frac[j][0] == 'Sediment/Detritus' and count == 0:
                    # TODO What really want here is weighted average of Ocs of regions maybe
                    total_nonlip_c += self.diet_frac[j][1] * (Ocs)
                    count += 1



        self.Vld = total_lip
        self.Vndc[region_index] = total_nonlip_c
        Vndm = total_nonlip_m
        self.Vwd[region_index] = 1 - (self.Vld + self.Vndc[region_index] + Vndm)


    #Where log is the dictonary of chemical concentrations for region
    def solve_steady_state(self, phi, i, Cwdp, Cwdo, log, chemical):

        denom = self.k_2[i] + self.k_e[i] + self.Kg\

        f_num = (self.k_1[i] * self.Mo * Cwdo) + (self.k_1[i] * self.Mp * Cwdp)

        l_num = 0

        for j in range (len(self.diet_frac)):
            if self.diet_frac[j][1] > 0:

                if self.diet_frac[j][0] == 'Sediment/Detritus':

                    concentration = chemical.Cs
                    l_num += (self.diet_frac[j][1]*concentration)

                else:
                    concentration = log[self.diet_frac[j][0]][chemical.name]
                    l_num += (self.diet_frac[j][1]*concentration)

        l_num = l_num * self.k_d[i]

        return (f_num + l_num)/denom


class Zooplank:


    def __init__(self, name, weight, vlb, flag, times, num_regions, num_chemicals, mp=0, vnb=.2, e_l=.72, e_n=.72, e_w=.25):

        self.name = name
        self.Wb = weight  # wet weight
        self.Vlb = vlb # Percent Lipid Content
        self.Vnb = vnb  # percent Nonlipid organic matter
        self.Mp = mp    # Percent Water Ventilated
        self.e_l = e_l  # Dietary absorption efficiency of lipids
        self.e_n = e_n  # Dietary absorption efficiency of nonlipid organic matter
        self.e_w = e_w  # Dietary absorption efficiency of water
        self.flag = flag
        self.Mo = 0  # Percent Over Water Ventilated
        self.Vwb = 0  # Percent Water Content
        self.Vld = 0  # lipid fraction of diet
        self.Vnd = 0  # Non-lipid fraction of diet
        self.Vwd = 0  # Water fraction of diet
        self.Vlg = 0  # lipid fraction of gut
        self.Vng = 0  # Non-lipid fraction of gut
        self.Vwg = 0  # water fraction of gut

        self.gd_set = 0
        self.kg_set = 0

        self.Gv = [[0 for _ in range(num_regions)] for _ in range (times)]    # Gill ventilation rate
        self.Gd = [[0 for _ in range(num_regions)] for _ in range (times)]    # Feeding rate
        self.Kg = [[0 for _ in range(num_regions)] for _ in range (times)]    # Growth rate
        self.Gf = [[0 for _ in range(num_regions)] for _ in range (times)]    # fecal egestion rate
        self.k_1 = [[[0 for _ in range (num_chemicals)] for _ in range(num_regions)] for _ in range (times)]    # clearance rate constant
        self.k_2 = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in range(times)]      # (rate constant chemical elem via respiratory)
        self.k_gb = [0 for _ in range (num_chemicals)]  # Gut biota partition coefficient
        self.k_e = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in range(times)]      # (rate constant via excretion into egested feces)
        self.k_d = [[[0 for _ in range(num_chemicals)] for _ in range(num_regions)] for _ in range(times)]      #(clearance rate constant via ingestion of food)

        self.Mb = 0  # mass of chemical for multiple times
        self.Cb = 0  # concentration for steady state

    def __str__(self):
        return str(self.__dict__)

    def set_el_en_ew(self, el, en, ew):
        self.e_l = el
        self.e_n = en
        self.e_w = ew

    def set_mp(self,mp):
        self.Mp = mp

    def set_gd(self,gd, times, len_reg):
        for i in range (times):
            for j in range(len_reg):
                self.Gd[i][j] = gd
        self.gd_set = 1

    def set_kg(self, kg, times, len_reg):

        for i in range (times):
            for j in range(len_reg):
                self.Kg[i][j] = kg
        self.kg_set = 1

    def set_vnb(self,vnb):
        self.Vnb = vnb

    def calc_vwb(self):
        self.Vwb = (1 - (self.Vlb + self.Vnb))

    def calc_mo(self):
        self.Mo = (1 - self.Mp)

    # sets the Gill ventilation rate of zooplank
    def calc_gv(self, regional_cox, time, region_index):
        try:
            self.Gv[time][region_index] = 1400 * ( math.pow(self.Wb, .65) / regional_cox[time])
        except ValueError:
            print('Weight of ' + self.name + ' was sampled as a negative number. This makes no sense. Check the range of the weight distribution. \n Caught at Gv calculation.')

    # sets G_d for this zooplank
    def calc_gd_filter(self, css, time_index ,region_index, sigma=1):
        self.Gd[time_index][region_index] = self.Gv[time_index][region_index] * css * sigma

    def calc_gd_no_filter(self, T, time_index, region_index):
        temp_at_time = T[time_index]
        self.Gd[time_index][region_index] = .022 * math.pow(self.Wb, .85) * math.exp((.06 * temp_at_time))

    def calc_kg(self, T_region, time_index, region_index):
        if T_region[time_index] < 15:
            self.Kg[time_index][region_index] = .0005 * math.pow(self.Wb, -.2)
        if T_region[time_index] >= 15:
            self.Kg[time_index][region_index] = .00251 * math.pow(self.Wb, -.2)

    def calc_gf(self, time_index, region_index):
        self.Gf[time_index][region_index] = (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) + ((1 - self.e_w) * self.Vwd)) * self.Gd[time_index][region_index]

    # returns k_1 of chemical for this zooplank
    def calc_k1(self, chem_ew, time_index, region_index, chemical_index):
        self.k_1[time_index][region_index][chemical_index] = chem_ew * (self.Gv[time_index][region_index] / self.Wb)

    # returns k_2 of chemical for this zooplank
    def calc_k2(self, chem_kow, k_1, time_index, region_index, chemical_index, beta=.035, density_lip=.9 , density_w=1):
        k_bw = (self.Vlb * chem_kow)/density_lip + (self.Vnb * beta * chem_kow) + (self.Vwb/density_w)
        self.k_2[time_index][region_index][chemical_index] = k_1[time_index][region_index][chemical_index] / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, ed, time_index, region_index, chemical_index):

        self.k_d[time_index][region_index][chemical_index] = ed * (self.Gd[time_index][region_index] / self.Wb)

    # sets the percentages of zooplank diet that are lipid, non-lipid and water
    def calc_diet_per(self, phyto):

        self.Vld = phyto.Vlb
        self.Vnd = phyto.Vnb
        self.Vwd = 1 - (self.Vld + self.Vnd)

    # sets the percentages of the gut
    def calc_gut_per(self):
        self.Vlg = ((1 - self.e_l) * self.Vld) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) + ((1 - self.e_w) * self.Vwd))
        self.Vng = ((1 - self.e_n) * self.Vnd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) + ((1 - self.e_w) * self.Vwd))
        self.Vwg = ((1 - self.e_w) * self.Vwd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) + ((1 - self.e_w) * self.Vwd))

    # returns K_gb for certian chemical
    def calc_kgb(self, chem_kow, chem_index,beta1=.35, beta2=.035, density_lip=.9, density_w=1, z_water=.05):
        top = ((self.Vlg * (z_water*chem_kow))/density_lip) + (self.Vng * beta1 * (z_water*chem_kow)) + (z_water*self.Vwg/density_w)
        bottom = ((self.Vlb * (z_water*chem_kow))/density_lip) + (self.Vnb * beta2 * z_water * chem_kow) + (z_water*self.Vwb/density_w)
        self.k_gb[chem_index] = top/bottom


    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, time_index, region_index,  chem_index):

        self.k_e[time_index][region_index][chem_index] = (self.k_gb[chem_index]/self.Wb) * chem_ed * self.Gf[time_index][region_index]


    # Where log is the dictonary of chemical concentrations
    def solve_steady_state(self, phi, i, Cwto, Cwdp,Cwdo, phyto_con):


        denom = self.k_2[i] + self.k_e[i] + self.Kg

        if phi != False:
            f_num = self.k_1[i] * (self.Mo * phi * Cwto + self.Mp * Cwdp)
        else:
            f_num = (self.k_1[i] * self.Mo * Cwdo) + (self.k_1[i] * self.Mp * Cwdp)
        f_num = f_num/1000

        l_num = phyto_con * self.k_d[i]

        return (f_num + l_num) / denom


class Pplank:

    # parameter #
    name = ''
    Vlb = None  # Percent Lipid Content
    Vnb = None  # Non-Lipid Organic Carbon Content
    Vwb = None  # Water Content
    Kg = None  # Growth-Rate
    A = None
    B = None

    # calc #
    k_1 = []  # List of k_1 (clearance rate constant) for each chemical. List should be number of chemicals long
    k_2 = []  # Same as k_1 but for k_2 (rate constant chemical elem via respiratory)

    # chemical
    Mb = None  # mass of chemical for multiple times
    Cb = None  # concentration for steady state

    def __init__(self, name, kg=None, a=.00006, b=5.5, vlp=.005, vnp=.065):
        self.name = name
        self.Vlb = vlp
        self.Vnb = vnp
        self.Kg = kg
        self.A = a
        self.B = b

    def __str__(self):
        return str(self.__dict__)

    def set_a_b(self,a,b):
        self.A = a
        self.B = b

    def set_vlp(self,vlp):
        self.Vlb = vlp

    def set_vnp(self,vnp):
        self.Vnb = vnp

    def calc_vwp(self):
        self.Vwb = 1 - (self.Vlb - self.Vnb)

    def set_kg(self, kg):
        self.Kg = kg


    # returns k_1 of chemical for this pplank dependent on chem kow
    def calc_k1(self, chem_kow):
        k1 = 1/(self.A + (self.B/chem_kow))
        return k1

    # returns k_2 of chemical for this pplank
    def calc_k2(self, chem_kow, k_1, density_lip=.9, density_w=1):
        k_pw = (self.Vlb * chem_kow) / density_lip + (self.Vnb * .35 * chem_kow) + (self.Vwb / density_w)
        return k_1 / k_pw


    def init_check(self):

        checks = 0

        # check 1
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry != None:
                count += 1

        if count == len(atts.values()):
            checks += 1

        #check 2
        if len(self.k_1) == len(self.k_2):
            checks += 1

        if checks == 2:
            return True

        return False

    def solve_steady_state(self, Cwd, i):

        return (self.k_1[i] * Cwd) / (self.k_2[i] + self.Kg)


class Chemical:

    # parameter #
    name = ''
    Kow = None  # Octanol water partition coefficient
    Cs = []  # Chemical concentration in sediment
    Cwto = []  # Chemical concentration in overlying water (total)
    Cwdo = []  # Chemical concentration in overlying water (dissolved)
    Ddoc = None  # Disequilibrium factor DOC
    Dpoc = None  # Disequilibrium factor POC

    # calculated #
    Ew = None  # Gill uptake effcifency
    Ed = None  # dietary chemical transfer efficiency
    phi = []  # Bioavailable solute fraction
    Cwp = []

    def __init__(self, name, kow, cs, len_regions):
        self.name = name
        self.Kow = math.pow(10,kow)
        self.Cs = cs
        self.phi = [0 for _ in range (len_regions)]
        self.Cwp = [0 for _ in range (len_regions)]
        # self.Cwto = cwto
        # self.Cwdo = cwdo
        # self.Ddoc = ddoc
        # self.Dpoc = dpoc
        self.calc_ew()
        self.calc_ed()

    def __str__(self):
        return self.__dict__.__str__()

    def set_cwto(self,cwto):
        self.Cwto = cwto

    def set_cwdo(self,cwdo):
        self.Cwdo = cwdo


    def set_ddoc_dpoc(self,ddoc,dpoc):
        self.Ddoc = ddoc
        self.Dpoc = dpoc


    def calc_ew(self):
        self.Ew = 1/(1.85 + (155/self.Kow))

    def calc_ed(self):
        self.Ed = math.pow((0.0000003*self.Kow + 2), -1)

    def calc_pore_water(self, region, index):
        Ocs = region.Ocs
        self.Cwp[index] = self.Cs[index]/(Ocs * .35 * self.Kow)

    def calc_phi_and_cwdo(self, region, index):

        if self.Cwto[index] != '' and self.Cwdo[index] != '':
            phi = self.Cwdo[index]/self.Cwto[index]
        if self.Cwto == '' and self.Cwdo != '':
            phi = False
        else:
            adoc = region.adoc
            apoc = region.apoc
            xdoc = region.Xdoc
            xpoc = region.Xpoc

            phi = 1/((1+xpoc*self.Dpoc*apoc*self.Kow)+(xdoc*self.Ddoc*adoc*self.Kow))


        self.phi[index] = phi

        if self.Cwdo[index] == ''and self.Cwto[index] != '':
            self.Cwdo[index] = (self.Cwto[index]/1000) * self.phi[index]

    def init_check(self):
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry != None:
                count += 1

        if count + 2 >= len(atts.values()):
            return True

        return False

class Region:

    def __init__(self, name, t, xdoc, xpoc, css, ocs, s, time_steps, adoc=.08, apoc=.35):

        self.name = name
        self.T = t   # mean water temperture
        self.Xdoc = xdoc  # Concentration Dissolved Organic Carbon Content
        self.Xpoc = xpoc  # Concentration Particulate Organic Carbon Content
        self.Css = css  # Concentration of suspended solids in water
        self.Ocs = ocs  # Sediment organic Carbon Content
        self.S = s  # Dissolved oxygen saturation
        self.Cox = [0 for _ in range (time_steps)]   # Dissolved Oxygen Concentration
        self.adoc = adoc  # DOC octanol proportionality constant
        self.apoc = apoc  # POC octanol proportionality constant

    def __str__(self):
        return str(self.__dict__)

    def set_adoc_apoc(self,apoc,adoc):
        self.adoc = adoc
        self.apoc = apoc

    def set_cox(self,cox, index):

            self.Cox[index] = cox

    def calc_cox(self, index):
        self.Cox[index] = (((-.24 * self.T[index]) + 14.04) * self.S)

    def init_check(self):
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry != None:
                count += 1

        if count == len(atts.values()):
            return True

        return False
