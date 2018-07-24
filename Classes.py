# The fish, invert, Pplank, Zooplank, and chemical objects has parameters for input that are from Arnot (2004) right now
import math


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
    Vnd = None  # Non-lipid fraction of diet
    Vwd = None  # Water fraction of diet
    Vlg = None  # lipid fraction of gut
    Vng = None  # Non-lipid fraction of gut
    Vwg = None  # water fraction of gut

    def __init__(self, name, weight, vlb, mp, diet_per, flag,vnb=.2, gd=None, e_l=.92, e_n=.6, e_w=.25, kg=None):
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


    # this could be improved...
    def set_vnb(self,vnb):
        self.Vnb = vnb
    def set_gd(self,gd):
        self.Gd = gd
    def set_el_en_ew(self,el,en,ew):
        self.e_l = el
        self.e_n = en
        self.e_w = ew

    def set_kg(self,kg):
        self.Kg = kg

    def calc_vwb(self):
        self.Vwb = (1 - (self.Vlb + self.Vnb))

    def calc_mo(self):
        self.Mo = (1 - self.Mp)

    ## Real calculations ##

    # sets the Gill ventilation rate of fish
    def calc_gv(self, region):
        region_cox = region.Cox
        self.Gv = 1400 * (math.pow(self.Wb, .65) / region_cox)

    # sets G_d for this fish
    def calc_gd(self, t, css,sigma=1):
        if self.flag == 1:
            self.Gd = self.Gv * css * sigma
        if self.flag == 0:
            self.Gd = .022 * math.pow(self.Wb, .85) * math.exp((.06 * t))

    def calc_kg(self, T):
        if T < 15:
            self.Kg = .0005 * math.pow(self.Wb, - .2)
        if T > 15:
            self.Kg = .00251 * math.pow(self.Wb, - .2)

    # returns k_1 of chemical for this fish
    def calc_k1(self, chem_ew):
        k1 = chem_ew * (self.Gv/self.Wb)
        return k1

    # returns k_2 of chemical for this fish
    def calc_k2(self, chem_kow, k_1, beta=.035):
        k_bw = (self.Vlb * chem_kow) + (self.Vnb * beta * chem_kow) + self.Vwb
        return k_1/k_bw

    # returns k_d of chemical for this fish
    def calc_kd(self, chem_ed):
        return chem_ed * (self.Gd/self.Wb)

    # sets the percentages of fish diet that are lipid, non-lipid and water
    def calc_diet_per(self, fishlog):

        total_nonlip = 0
        total_lip = 0
        total_all = 0

        # loop through possible types of fishes
        for i in range(len(fishlog)):
            # loop through fish names in diet
            for j in range(len(self.diet_frac)):
                # finding right fish
                if fishlog[i].name == self.diet_frac[j][0]:


                    l_toadd = self.diet_frac[j][1] * (fishlog[i].Vlb * fishlog[i].Wb)
                    nl_toadd = self.diet_frac[j][1] * (fishlog[i].Vnb * fishlog[i].Wb)

                    total_nonlip += nl_toadd
                    total_lip += l_toadd
                    total_all += (self.diet_frac[j][1] * fishlog[i].Wb)

        self.Vld = total_lip/total_all
        self.Vnd = total_nonlip/total_all
        self.Vwd = 1 - (self.Vld + self.Vnd)

    # sets the percentages of the gut
    def calc_gut_per(self):
        self.Vlg = ((1-self.e_l)*self.Vld)/(((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vnd) + ((1-self.e_w)*self.Vwd))
        self.Vng = ((1-self.e_n)*self.Vnd)/(((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vnd) + ((1-self.e_w)*self.Vwd))
        self.Vwg = ((1-self.e_w)*self.Vwd)/(((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vnd) + ((1-self.e_w)*self.Vwd))

    # sets Gf for this zooplank
    def calc_gf(self):
        self.Gf = (((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vnd)+((1-self.e_w)*self.Vwd))*self.Gd

    # returns K_gb for certian chemical
    def calc_kgb(self, chem_kow, beta=.035):
        k_gb = (self.Vlg*chem_kow+self.Vng*beta*chem_kow+self.Vng)/(self.Vlb*chem_kow+self.Vnb*beta*chem_kow+self.Vwb)
        return k_gb

    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, k_gb):
        k_e = self.Gf * chem_ed * (k_gb/self.Wb)
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
    def solve_steady_state(self, phi, i, Cwto, Cwds, log, chemical):

        denom = self.k_2[i] + self.k_e[i] + self.Kg
        f_num = self.k_1[i] * (self.Mo * phi * Cwto + self.Mp * Cwds)

        l_num = 0

        for j in range (len(self.diet_frac)):
            if self.diet_frac[j][1] > 0:
                concentration = log[self.diet_frac[j][0]][chemical.name]
                l_num += (self.diet_frac[j][1]*concentration)

        l_num = l_num * self.k_d[i]

        return (f_num+l_num)/denom


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

    def __init__(self, name, weight, lip_con, mp=0, vnb=.2, gd=None, e_l=.72, e_n=.72, e_w=.25, kg=None):

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
        self.Gv = 1400 * (math.pow(self.Wb, .65) / region_cox)

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
    def calc_k2(self, chem_kow, k_1, beta=.035):
        k_bw = (self.Vlb * chem_kow) + (self.Vnb * beta * chem_kow) + self.Vwb
        return k_1 / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, chem_ed):
        return chem_ed * (self.Gd / self.Wb)

    # sets the percentages of zooplank diet that are lipid, non-lipid and water
    def calc_diet_per(self, phyto):

        self.Vld = phyto.Vlp
        self.Vnd = phyto.Vnp
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
    def calc_kgb(self, chem_kow, beta=.035):
        k_gb = (self.Vlg * chem_kow + self.Vng * beta * chem_kow + self.Vng) / (self.Vlb * chem_kow + self.Vnb * beta * chem_kow + self.Vwb)
        return k_gb

    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, k_gb):
        k_e = self.Gf * chem_ed * (k_gb / self.Wb)
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

        l_num = phyto_con * self.k_d[i]

        return (f_num + l_num) / denom


class Pplank:

    # parameter #
    name = ''
    Vlp = None  # Percent Lipid Content
    Vnp = None  # Non-Lipid Organic Carbon Content
    Vwp = None  # Water Content
    Kg = None  # Growth-Rate
    A = None
    B = None

    # calc #
    k_1 = []  # List of k_1 (clearance rate constant) for each chemical. List should be number of chemicals long
    k_2 = []  # Same as k_1 but for k_2 (rate constant chemical elem via respiratory)

    def __init__(self, name, kg, a=.00006, b=5.5, vlp=.005, vnp=.065):
        self.name = name
        self.Vlp = vlp
        self.Vnp = vnp
        self.Kg = kg
        self.A = a
        self.B = b


    def set_a_b(self,a,b):
        self.A = a
        self.B = b

    def set_vlp(self,vlp):
        self.Vlp = vlp

    def set_vnp(self,vnp):
        self.Vnp = vnp

    def calc_vwp(self):
        self.Vwp = 1 - (self.Vlp - self.Vnp)

    # returns k_1 of chemical for this pplank dependent on chem kow
    def calc_k1(self, chem_kow):
        k1 = math.pow((self.A + self.B)/chem_kow, -1)
        return k1

    # returns k_2 of chemical for this pplank
    def calc_k2(self, chem_kow, k_1):
        k_pw = (self.Vlp * chem_kow) + (self.Vnp * .35 * chem_kow) + self.Vwp
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
    Ddoc = None  # Disequilibrium factor DOC
    Dpoc = None  # Disequilibrium factor POC

    # calculated #
    Ew = None  # Gill uptake effcifency
    Ed = None  # dietary chemical transfer efficiency
    phi = None  # Bioavailable solute fraction

    def __init__(self, name, kow, cs, cwto=None, cwdo=None, ddoc=None, dpoc=None):
        self.name = name
        self.Kow = kow
        self.Cs = cs
        self.Cwto = cwto
        self.Cwdo = cwdo
        self.Ddoc = ddoc
        self.Dpoc = dpoc
        self.calc_ew()
        self.calc_ed()

    def set_cwto(self,cwto):
        self.Cwto = cwto

    def set_cwdo(self,cwdo):
        self.Cwdo = cwdo


    def set_ddoc_dpoc(self,ddoc,dpoc):
        self.Ddoc = ddoc
        self.Dpoc = dpoc


    def calc_ew(self):
        self.Ew = (1.85 + math.pow((155/self.Kow), -1))

    def calc_ed(self):
        self.Ed = math.pow((0.0000003*self.Kow + 2), -1)

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

            self.Cwdo = self.Cwto * self.phi

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

    def set_adoc_apoc(self,apoc,adoc):
        self.adoc = adoc
        self.apoc = apoc

    def set_cox(self,cox):
        self.Cox = cox

    def calc_cox(self):
        self.Cox = (-.024 * self.T + 14.04) * self.S

    def init_check(self):
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry != None:
                count += 1

        if count == len(atts.values()):
            return True

        return False


