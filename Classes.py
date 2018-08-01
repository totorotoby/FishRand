# The fish, invert, Pplank, Zooplank, and chemical objects has parameters for input that are from Arnot (2004) right now
import math
#
class Var:

    def __init__(self, type, dist, param):

        self.type = type   # 0 is uncertain , 1 is variable, 2 is both, 3 is point estimate
        self.dist = dist   # distribution name
        self.param = param
        self.value = None
        self.lhs = None

    def __str__(self):

        to_print = '\ntype | ' + self.type + '\ndistribution | ' + self.dist +'\nparameters | ' + str(self.param)
        return to_print

class Fish:
    # parameter #
    name = ''
    Wb = None  # wet weight
    Vlb = None  # percent Lipid Content
    Vnb = None  # percent Nonlipid organic matter
    Mp = None  # Percent Pore Water Ventilated
    diet_frac = None # List frac of diet. The sum of entries must add up to 1, and Length should be number of organisms
    e_l = None  # Dietary absorption efficiency of lipid
    e_n = None  # Dietary absorption efficiency of nonlipid organic matter
    e_w = None  # Dietary absorption efficiency of water
    flag = 0  # Filter feeder flag

    # parameter or calc #
    Gd = None  # Feeding rate
    Kg = None  # Growth rate (only for steady state)

    # calc #
    Vwb = None  # Percent of water content (1 - (Vlb + Vnb))
    Mo = None  # Percent Over Water Ventilated
    k_1 = None  # List of k_1 (clearance rate constant via resp) for each chemical.List should be number of chemicals long
    k_2 = None  # Same as k_1 but for k_2 (rate constant chemical elem via respiratory)
    k_d = None  # Same but for k_d (clearance rate constant via ingestion of food)
    k_e = None  # (rate constant via excretion into egested feces)
    k_gb = None  # Gut–biota partition coefficient
    Gv = None  # Gill ventilation rate
    Gf = None  # fecal egestion rate
    Vld = None  # lipid fraction of diet
    Vndc = None  # Non-lipid fraction of diet
    Vndm = None
    Vwd = None  # Water fraction of diet
    Vlg = None  # lipid fraction of gut
    Vngc = None  # Non-lipid fraction of gut
    Vngm = None
    Vwg = None  # water fraction of gut

    # chemical
    Mb = None  # mass of chemical for multiple times
    Cb = None  # concentration for steady state


    def __init__(self, name, weight, vlb, mp, diet_per, flag,vnb=.2, gd=None, e_l=.75, e_n=.75, e_w=.5, kg=None):
        self.name = name
        self.Wb = weight
        self.Vlb = vlb
        self.Vnb = vnb
        self.Mp = mp
        self.Gd = gd
        self.diet_frac = diet_per
        self.e_l = e_l
        self.e_n = e_n
        self.e_w = e_w
        self.Kg = kg
        self.flag = flag
        self.Cb = []

    def __str__(self):
        return self.__dict__.__str__()

    def set_el_en_ew(self, el, en, ew):
        self.e_l = el
        self.e_n = en
        self.e_w = ew

    def set_mp(self,mp):
        self.Mp = mp

    def set_gd(self,gd):
        self.Gd = gd

    def set_kg(self, kg):
        self.Kg = kg

    def set_vnb(self,vnb):
        self.Vnb = vnb

    def calc_vwb(self):
        self.Vwb = (1 - (self.Vlb + self.Vnb))

    def calc_mo(self):
        self.Mo = (1 - self.Mp)

    # sets the Gill ventilation rate of zooplank
    # TODO there was a domain error with one of the GV's has to do with negatives?
    def calc_gv(self, region):
        region_cox = region.Cox
        self.Gv = 1400 * ( math.pow(self.Wb, .65) / region_cox)

    # sets G_d for this zooplank
    def calc_gd(self, t, css, flag,sigma=1):
        if flag == 1:
            self.Gd = self.Gv * css * sigma
        if flag == 0:
            self.Gd = .022 * math.pow(self.Wb, .85) * math.exp((.06 * t))

    def calc_kg(self, T):
        if T < 15:
            self.Kg = .0005 * math.pow(self.Wb, -.2)
        if T > 15:
            self.Kg = .00251 * math.pow(self.Wb, -.2)

    # returns k_1 of chemical for this zooplank
    def calc_k1(self, chem_ew):
        k1 = chem_ew * (self.Gv / self.Wb)
        return k1

    # returns k_2 of chemical for this zooplank
    def calc_k2(self, chem_kow, k_1, beta=.035, density_lip=.9 , density_w=1):

        k_bw = (self.Vlb * chem_kow)/density_lip + (self.Vnb * beta * chem_kow) + (self.Vwb/density_w)
        return k_1 / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, chem_ed):

        return chem_ed * (self.Gd / self.Wb)

    # sets the percentages of zooplank diet that are lipid, non-lipid and water
    def calc_diet_per(self, fishlog, region):

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

                    total_nonlip_c += self.diet_frac[j][1] * (region.Ocs)
                    count += 1



        self.Vld = total_lip
        self.Vndc = total_nonlip_c
        self.Vndm = total_nonlip_m
        self.Vwd = 1 - (self.Vld + self.Vndc + self.Vndm)


    # sets the percentages of the gut
    def calc_gut_per(self):

        self.Vlg = ((1 - self.e_l) * self.Vld) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vngc = ((1 - self.e_n) * self.Vndc) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vngm = ((1 - self.e_n) * self.Vndm) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vwg = ((1 - self.e_w) * self.Vwd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc)+ ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))

    # sets Gf for this zooplank
    def calc_gf(self):

        self.Gf = (((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vndc) + ((1-self.e_n)*self.Vndm) + ((1-self.e_w)*self.Vwd))*self.Gd

    # returns K_gb for certian chemical
    def calc_kgb(self, chem_kow, beta1=.35, beta2=.035, density_lip=.9, density_w=1, z_water=.05):

        top = ((self.Vlg * (z_water*chem_kow))/density_lip) + (self.Vngc * beta1 * (z_water*chem_kow)) + (self.Vngm * beta2 * (z_water*chem_kow)) + (z_water*self.Vwg/density_w)
        bottom = ((self.Vlb * (z_water*chem_kow))/density_lip) + (self.Vnb * beta2 * z_water * chem_kow) + (z_water*self.Vwb/density_w)
        k_gb = top/bottom
        return k_gb

    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, k_gb):
        k_e = (k_gb/self.Wb) * chem_ed * self.Gf
        return k_e

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

    #Where log is the dictonary of chemical concentrations for region
    def solve_steady_state(self, phi, i, Cwp, Cwds, log, chemical):

        denom = self.k_2[i] + self.k_e[i] + self.Kg
        f_num = (self.k_1[i] * self.Mo * Cwds) + (self.k_1[i] * self.Mp * Cwp)

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

    # parameter #
    name = ''
    Wb = None  # wet weight
    Vlb = None  # Percent Lipid Content
    Vnb = None  # percent Nonlipid organic matter
    Mp = None  # Percent Pore Water Ventilated
    diet_frac = None  # List of fractions of diet. The sum of entries must add up to 1
    e_l = None  # Dietary absorption efficiency of lipids
    e_n = None  # Dietary absorption efficiency of nonlipid organic matter
    e_w = None  # Dietary absorption efficiency of water
    flag = 0

    # parameter or calc #
    Gd = None  # Feeding rate
    Kg = None  # Growth rate

    # calc #
    Mo = None  # Percent Over Water Ventilated
    Vwb = None  # Percent Water Content
    k_1 = []  # List of k_1 (clearance rate constant) for each chemical. List should be number of chemicals long
    k_2 = []  # Same as k_1 but for k_2 (rate constant chemical elem via respiratory)
    k_d = []  # Same but for k_d (clearance rate constant via ingestion of food)
    k_e = []  # (rate constant via excretion into egested feces)
    k_gb = []  # Gut–biota partition coefficient
    Gv = None  # Gill ventilation rate
    Gf = None  # fecal egestion rate
    Vld = None  # lipid fraction of diet
    Vnd = None  # Non-lipid fraction of diet
    Vwd = None  # Water fraction of diet
    Vlg = None  # lipid fraction of gut
    Vng = None  # Non-lipid fraction of gut
    Vwg = None  # water fraction of gut

    # chemical
    Mb = None  # mass of chemical for multiple times
    Cb = None  # concentration for steady state

    def __init__(self, name, weight, lip_con, flag, mp=0, vnb=.2, gd=None, e_l=.72, e_n=.72, e_w=.25, kg=None):

        self.name = name
        self.Wb = weight
        self.Vlb = lip_con
        self.Vnb = vnb
        self.Mp = mp
        self.Gd = gd
        self.Kg = kg
        self.e_l = e_l
        self.e_n = e_n
        self.e_w = e_w
        self.flag = flag

    def __str__(self):
        return self.__dict__

    def set_el_en_ew(self, el, en, ew):
        self.e_l = el
        self.e_n = en
        self.e_w = ew

    def set_mp(self,mp):
        self.Mp = mp

    def set_gd(self,gd):
        self.Gd = gd

    def set_kg(self, kg):
        self.Kg = kg

    def set_vnb(self,vnb):
        self.Vnb = vnb

    def calc_vwb(self):
        self.Vwb = (1 - (self.Vlb + self.Vnb))

    def calc_mo(self):
        self.Mo = (1 - self.Mp)

    # sets the Gill ventilation rate of zooplank
    def calc_gv(self, region):
        region_cox = region.Cox
        self.Gv = 1400 * ( math.pow(self.Wb, .65) / region_cox)

    # sets G_d for this zooplank
    def calc_gd(self, t, css, flag,sigma=1):
        if flag == 1:
            self.Gd = self.Gv * css * sigma
        if flag == 0:
            self.Gd = .022 * math.pow(self.Wb, .85) * math.exp((.06 * t))

    def calc_kg(self, T):
        if T < 15:
            self.Kg = .0005 * math.pow(self.Wb, -.2)
        if T > 15:
            self.Kg = .00251 * math.pow(self.Wb, -.2)

    # returns k_1 of chemical for this zooplank
    def calc_k1(self, chem_ew):
        k1 = chem_ew * (self.Gv / self.Wb)
        return k1

    # returns k_2 of chemical for this zooplank
    def calc_k2(self, chem_kow, k_1, beta=.035, density_lip=.9 , density_w=1):
        k_bw = (self.Vlb * chem_kow)/density_lip + (self.Vnb * beta * chem_kow) + (self.Vwb/density_w)
        return k_1 / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, chem_ed):

        return chem_ed * (self.Gd / self.Wb)

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
    # sets Gf for this zooplank
    def calc_gf(self):
        self.Gf = (((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vnd)+((1-self.e_w)*self.Vwd))*self.Gd

    # returns K_gb for certian chemical
    def calc_kgb(self, chem_kow, beta1=.35, beta2=.035, density_lip=.9, density_w=1, z_water=.05):
        top = ((self.Vlg * (z_water*chem_kow))/density_lip) + (self.Vng * beta1 * (z_water*chem_kow)) + (z_water*self.Vwg/density_w)
        bottom = ((self.Vlb * (z_water*chem_kow))/density_lip) + (self.Vnb * beta2 * z_water * chem_kow) + (z_water*self.Vwb/density_w)
        k_gb = top/bottom
        return k_gb

    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, k_gb):

        k_e = (k_gb/self.Wb) * chem_ed * self.Gf
        return k_e

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

    # Where log is the dictonary of chemical concentrations
    def solve_steady_state(self, phi, i, Cwto, Cwds, phyto_con):


        denom = self.k_2[i] + self.k_e[i] + self.Kg

        f_num = self.k_1[i] * (self.Mo * phi * Cwto + self.Mp * Cwds)
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

    def __init__(self, name, kg, a=.00006, b=5.5, vlp=.005, vnp=.065):
        self.name = name
        self.Vlb = vlp
        self.Vnb = vnp
        self.Kg = kg
        self.A = a
        self.B = b

    def __str__(self):
        return self.__dict__.__str__()

    def set_a_b(self,a,b):
        self.A = a
        self.B = b

    def set_vlp(self,vlp):
        self.Vlb = vlp

    def set_vnp(self,vnp):
        self.Vnb = vnp

    def calc_vwp(self):
        self.Vwb = 1 - (self.Vlb - self.Vnb)

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
    Kow = None  # Octanol–water partition coefficient
    Cs = None  # Chemical concentration in sediment
    Cwto = None  # Chemical concentration in overlying water (total)
    Cwdo = None  # Chemical concentration in overlying water (dissolved)
    Cwp = None
    Ddoc = None  # Disequilibrium factor DOC
    Dpoc = None  # Disequilibrium factor POC

    # calculated #
    Ew = None  # Gill uptake effcifency
    Ed = None  # dietary chemical transfer efficiency
    phi = None  # Bioavailable solute fraction

    def __init__(self, name, kow, cs, cwto=None, cwdo=None, ddoc=None, dpoc=None):
        self.name = name
        self.Kow = math.pow(10,kow)
        self.Cs = cs
        self.Cwto = cwto
        self.Cwdo = cwdo
        self.Ddoc = ddoc
        self.Dpoc = dpoc
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

    def calc_pore_water(self, region):
        Ocs = region.Ocs
        self.Cwp = self.Cs/(Ocs * .35 * self.Kow)


    def calc_phi_and_cwdo(self, region):

        if self.Cwto and self.Cwdo != None:
            self.phi = self.Cwdo/self.Cwto
        else:
            adoc = region.adoc
            apoc = region.apoc
            xdoc = region.Xdoc
            xpoc = region.Xpoc

            self.phi = 1/((1+xpoc*self.Dpoc*apoc*self.Kow)+(xdoc*self.Ddoc*adoc*self.Kow))

        if self.Cwdo == None:

            self.Cwdo = (self.Cwto/1000) * self.phi

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

    # parameter #
    name = ''
    T = None  # mean water temperture
    Xdoc = None  # Concentration Dissolved Organic Carbon Content
    Xpoc = None  # Concentration Particulate Organic Carbon Content
    Css = None  # Concentration of suspended solids in water
    Ocs = None  # Sediment organic Carbon Content
    S = None  # Dissolved oxygen saturation
    adoc = None  # DOC–octanol proportionality constant
    apoc = None  # POC–octanol proportionality constant

    # parameter or calc #
    Cox = None  # Dissolved Oxygen Content

    def __init__(self, name, t, xdoc, xpoc, css, ocs, s, cox=None, adoc=.08, apoc=.35):

        self.T = t
        self.name = name
        self.Xdoc = xdoc
        self.Xpoc = xpoc
        self.Css = css
        self.Ocs = ocs
        self.S = s
        self.Cox = cox
        self.adoc = adoc
        self.apoc = apoc

    def __str__(self):
        return self.__dict__

    def set_adoc_apoc(self,apoc,adoc):
        self.adoc = adoc
        self.apoc = apoc

    def set_cox(self,cox):
        self.Cox = cox

    def calc_cox(self):
        self.Cox = ((-.24 * self.T) + 14.04) * self.S

    def init_check(self):
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry != None:
                count += 1

        if count == len(atts.values()):
            return True

        return False


