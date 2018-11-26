import math
import Bioaccum

class Zooplank:

    def __init__(self, name, weight, vlb, flag, num_chemicals, per_step, vnb=.2, e_l=.72, e_n=.72, e_w=.25):
        self.name = name
        self.Wb = weight  # wet weight
        self.Vlb = vlb  # Percent Lipid Content
        self.Vnb = vnb  # percent Nonlipid organic matter
        self.Mp = 0  # Percent Water Ventilated
        self.e_l = e_l  # Dietary absorption efficiency of lipids
        self.e_n = e_n  # Dietary absorption efficiency of nonlipid organic matter
        self.e_w = e_w  # Dietary absorption efficiency of water
        self.flag = flag
        self.Mo = 0  # Percent Over Water Ventilated
        self.Vwb = 0  # Percent Water Content
        self.Vld = 0  # lipid fraction of diet
        self.Vwd = 0  # Water fraction of diet
        self.Vlg = 0  # lipid fraction of gut
        self.Vwg = 0  # water fraction of gut
        self.Vndc = 0  # nonlip fraction of diet from dirt
        self.Vndm = 0  # nonlip fracion of diet from animal
        self.Vngc = 0  # Non-lipid fraction of gut
        self.Vngm = 0

        self.gd_set = 0
        self.kg_set = 0

        self.Gv = 0  # Gill ventilation rate
        self.Gd = 0  # Feeding rate
        self.Kg = 0  # Growth rate
        self.Gf = 0  # fecal egestion rate
        self.k_1 = [0 for _ in range(num_chemicals)]  # clearance rate constant
        self.k_2 = [0 for _ in range(num_chemicals)]  # (rate constant chemical elem via respiratory)
        self.k_gb = [0 for _ in range(num_chemicals)]  # Gut biota partition coefficient
        self.k_e = [0 for _ in range(num_chemicals)]  # (rate constant via excretion into egested feces)
        self.k_d = [0 for _ in range(num_chemicals)]  # (clearance rate constant via ingestion of food)

        self.days_per_step = per_step

    # a display but only for the first chemical
    def __str__(self):

        string = '\n'
        string += str(self.name)
        string += '\n'
        string += '|         K_1         |         K_D         |         K_2         |         K_E       |         K_GB       |         K_g       |'
        string += '\n'
        string += '|  ' + str(round(self.k_1[0],16)) + '  |'
        string += '  ' + str(round(self.k_d[0],16)) + '  |'
        string += '  ' + str(round(self.k_2[0],16)) + '  |'
        string += '  ' + str(round(self.k_e[0],16)) + '  |'
        string += '  ' + str(round(self.k_gb[0], 16)) + '  |'
        string += '  ' + str(round(self.Kg, 16)) + '  |'
        string += '\n\n'
        string += 'm_o :  ' + str(self.Mo) + '\n'
        string += 'm_p :  ' + str(self.Mp) + '\n\n\n'


        return string

    def set_el_en_ew(self, el, en, ew):
        self.e_l = el
        self.e_n = en
        self.e_w = ew

    def set_mp(self, mp):
        self.Mp = mp

    def set_gd(self, gd):
        self.Gd = gd

    def set_kg(self, kg):
        self.Kg = kg

    def set_vnb(self, vnb):
        self.Vnb = vnb

    def calc_vwb(self):
        self.Vwb = (1 - (self.Vlb + self.Vnb))

    def calc_mo(self):
        self.Mo = (1 - self.Mp)

    # sets the Gill ventilation rate of zooplank
    def calc_gv(self, region_cox):
        self.Gv = 1400 * (math.pow(self.Wb, .65) / region_cox)

    # sets G_d for this zooplank
    def calc_gd_filter(self, css, sigma=1):
        
        self.Gd = self.Gv * css * sigma

    def calc_gd_no_filter(self, T):

        self.Gd = .022 * math.pow(self.Wb, .85) * math.exp((.06 * T))

    def calc_kg(self, T):

        if T <= 18:
            self.Kg = .000502 * math.pow(self.Wb, -.2)
        if T > 18:
            self.Kg = .00251 * math.pow(self.Wb, -.2)

    # returns k_1 of chemical for this zooplank
    def calc_k1(self, chem_ew, chem_index):
        self.k_1[chem_index] = chem_ew * (self.Gv / self.Wb)

    # returns k_2 of chemical for this zooplank
    def calc_k2(self, chem_kow, chem_index, beta1, beta5, density_lip=.9, density_w=1):




        k_bw = (self.Vlb * beta1 * chem_kow)/density_lip + (self.Vnb * beta5 * chem_kow) + (self.Vwb/density_w)



        self.k_2[chem_index] = self.k_1[chem_index] / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, chem_ed, chem_index):

        self.k_d[chem_index] = chem_ed * (self.Gd / self.Wb)

    # sets the percentages of zooplank diet that are lipid, non-lipid and water
    def calc_diet_per(self, phyto):

        self.Vld = phyto.Vlb
        self.Vnd = phyto.Vnb
        self.Vwd = 1 - (self.Vld + self.Vnd)

    # sets the percentages of the gut
    def calc_gut_per(self):
        self.Vlg = ((1 - self.e_l) * self.Vld) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) +
                                                  ((1 - self.e_w) * self.Vwd))
        self.Vng = ((1 - self.e_n) * self.Vnd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) +
                                                  ((1 - self.e_w) * self.Vwd))
        self.Vwg = ((1 - self.e_w) * self.Vwd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vnd) +
                                                  ((1 - self.e_w) * self.Vwd))
    # sets Gf for this zooplank

    def calc_gf(self):

        self.Gf = (((1-self.e_l)*self.Vld) + ((1-self.e_n)*self.Vnd)+((1-self.e_w)*self.Vwd))*self.Gd

    # returns K_gb for certian chemical

    def calc_kgb(self, chem_kow, chem_index, beta3, beta4 ,beta5, density_lip=.9, density_w=1, z_water=.05):


        top = ((self.Vlg * (z_water*beta3*chem_kow))/density_lip) + (self.Vng * beta4 * (z_water*chem_kow)) +\
              (z_water*self.Vwg/density_w)

        bottom = ((self.Vlb * (z_water*beta3*chem_kow))/density_lip) + (self.Vnb * beta5 * z_water * chem_kow) +\
                 (z_water*self.Vwb/density_w)

        k_gb = top/bottom


        self.k_gb[chem_index] = k_gb


    # returns k_e for certain chemical
    def calc_ke(self, chem_ed, chem_index):
        k_e = (self.k_gb[chem_index]/self.Wb) * chem_ed * self.Gf
        self.k_e[chem_index] = k_e


    def init_check(self):

        checks = 0

        # check 1
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry is not None:
                count += 1

        if count == len(atts.values()):
            checks += 1

        # check 2
        if len(self.k_1) == len(self.k_2):
            checks += 1

        if checks == 2:
            return True

        return False

    # Where log is the dictonary of chemical concentrations
    def solve_steady_state(self, phi, chem_index, Cwp, Cwdo, phyto_con):

        denom = self.k_2[chem_index] + self.k_e[chem_index] + self.Kg

        f_num = (self.k_1[chem_index] * self.Mo * Cwdo) + (self.k_1[chem_index] * self.Mp * Cwp)


        l_num = phyto_con * self.k_d[chem_index]

        return (f_num + l_num) / denom


    def solve_next_time_step(self, phi, chem_index, Cwp, Cwdo, phyto_con, pre_step):

        f_num = (self.k_1[chem_index] * self.Mo * Cwdo) + (self.k_1[chem_index] * self.Mp * Cwp)

        l_num = phyto_con * self.k_d[chem_index]

        q = f_num + l_num
        k = self.k_2[chem_index] + self.k_e[chem_index]

        top = (pre_step * ((1 / self.days_per_step) - (k / 2)) + q)
        bottom = (1 / self.days_per_step) + (k / 2)

        return top/bottom


class Fish(Zooplank):

    def __init__(self, name, weight, vlb, diet_data, flag, num_chemicals, per_step, vnb=.2, e_l=.75, e_n=.75, e_w=.5):
        Zooplank.__init__(self, name, weight, vlb, flag, num_chemicals, per_step, vnb, e_l, e_n, e_w)
        self.diet_frac = diet_data

    # sets the percentages of zooplank diet that are lipid, non-lipid and water
    def calc_diet_per(self, fishlog, region=None):

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

                    l_toadd = self.diet_frac[j][1] * fishlog[i].Vlb
                    nl_toadd = self.diet_frac[j][1] * fishlog[i].Vnb

                    total_nonlip_c += nl_toadd
                    total_lip += l_toadd

                elif fishlog[i].name == self.diet_frac[j][0]:

                    l_toadd = self.diet_frac[j][1] * fishlog[i].Vlb
                    nl_toadd = self.diet_frac[j][1] * fishlog[i].Vnb

                    total_nonlip_m += nl_toadd
                    total_lip += l_toadd

                elif self.diet_frac[j][0] == 'Sediment/Detritus' and self.diet_frac[j][1] > 0  and count == 0:

                    total_nonlip_c += self.diet_frac[j][1] * region.Ocs
                    count += 1

        self.Vld = total_lip
        self.Vndc = total_nonlip_c
        self.Vndm = total_nonlip_m
        self.Vwd = 1 - (self.Vld + self.Vndc + self.Vndm)

    # sets the percentages of the gut
    def calc_gut_per(self):

        self.Vlg = ((1 - self.e_l) * self.Vld) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc) +
                                                  ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vngc = ((1 - self.e_n) * self.Vndc) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc) +
                                                    ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vngm = ((1 - self.e_n) * self.Vndm) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc) +
                                                    ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))
        self.Vwg = ((1 - self.e_w) * self.Vwd) / (((1 - self.e_l) * self.Vld) + ((1 - self.e_n) * self.Vndc) +
                                                  ((1 - self.e_n) * self.Vndm) + ((1 - self.e_w) * self.Vwd))

    # sets Gf for this zooplank
    def calc_gf(self):

        self.Gf = (((1-self.e_l) * self.Vld) + ((1-self.e_n) * self.Vndc) + ((1-self.e_n) * self.Vndm) +
                   ((1-self.e_w) * self.Vwd)) * self.Gd

    # returns K_gb for certian chemical
    def calc_kgb(self, chem_kow, chem_index, beta3, beta4, beta5, density_lip=.9, density_w=1, z_water=.05):

        top = ((self.Vlg * (z_water*beta3*chem_kow))/density_lip) + (self.Vngc * beta4 * (z_water*chem_kow)) +\
              (self.Vngm * beta5 * (z_water*chem_kow)) + (z_water*self.Vwg/density_w)
        bottom = ((self.Vlb * (z_water*beta3*chem_kow))/density_lip) + (self.Vnb * beta5 * z_water * chem_kow) +\
                 (z_water*self.Vwb/density_w)
        k_gb = top/bottom

        self.k_gb[chem_index] = k_gb


    #Where log is the dictonary of chemical concentrations for region
    def solve_steady_state(self, phi, chem_index, Cwp, Cwdo, log, chemical):


        denom = self.k_2[chem_index] + self.k_e[chem_index] + self.Kg
        f_num = (self.k_1[chem_index] * self.Mo * Cwdo) + (self.k_1[chem_index] * self.Mp * Cwp)


        l_num = 0

        for j in range(len(self.diet_frac)):
            if self.diet_frac[j][1] > 0:

                if self.diet_frac[j][0] == 'Sediment/Detritus':

                    concentration = chemical.Cs
                    l_num += (self.diet_frac[j][1]*concentration)

                else:
                    concentration = log[self.diet_frac[j][0]][chemical.name]
                    l_num += (self.diet_frac[j][1]*concentration)

       

        l_num = l_num * self.k_d[chem_index]

        return (f_num + l_num)/denom

    def solve_next_time_step(self, phi, chem_index, Cwp, Cwdo, fishlog_new, fishlog_old, chemical, pre_step, out_check):

        if out_check == 0:


            k = self.k_2[chem_index] + self.k_e[chem_index]

            q1 = (self.k_1[chem_index] * self.Mo * Cwdo) + (self.k_1[chem_index] * self.Mp * Cwp)

            q2 = 0

            if Bioaccum.not_eating(self) == True:

                q2 = 0


            else:
                for j in range(len(self.diet_frac)):

                    if self.diet_frac[j][1] > 0:

                        if self.diet_frac[j][0] == 'Sediment/Detritus':

                            concentration = chemical.Cs
                            q2 += (self.diet_frac[j][1]*concentration)

                        else:
                            if type(fishlog_old) == list:

                                try:
                                    old_conc = fishlog_old[0][self.diet_frac[j][0]][chemical.name]
                                    new_conc = fishlog_new[0][self.diet_frac[j][0]][chemical.name]

                                except:
                                    if fishlog_old[1] == None:
                                        old_conc = 0
                                    else:
                                        old_conc = fishlog_old[1][self.diet_frac[j][0]][chemical.name]
                                    new_conc = fishlog_new[1][self.diet_frac[j][0]][chemical.name]
                            else:
                                old_conc = fishlog_old[self.diet_frac[j][0]][chemical.name]
                                new_conc = fishlog_new[self.diet_frac[j][0]][chemical.name]

                            if old_conc == 0:
                                concentration = new_conc
                            else:
                                concentration = (old_conc + new_conc)/2
                            q2 += (self.diet_frac[j][1]*concentration)

            q2 = q2 * self.k_d[chem_index]


            q = q1 + q2

            top = (pre_step * ((1 / self.days_per_step) - (k / 2)) + q)
            bottom = (1 / self.days_per_step) + (k / 2)

            change = ((top/bottom) - pre_step)

        # outside of pond so only take the subtraction part
        else:
            if pre_step == 0:
                change = 0
            else:
                change = (self.k_2[chem_index] + self.k_e[chem_index]) * pre_step

        return change


class Pplank:

    def __init__(self, name, kg, len_chem, per_step, a=.00006, b=5.5, vlb=.005, vnb=.065):
        self.name = name
        self.Vlb = vlb
        self.Vnb = vnb
        self.Vwb = 0
        self.Kg = kg
        self.A = a
        self.B = b
        self.k_1 = [0 for _ in range(len_chem)]  # List of k_1 (clearance rate constant) for each chemical.
        self.k_2 = [0 for _ in range(len_chem)]  # Same as k_1 but for k_2
        #  (rate constant chemical elem via respiratory)

        self.days_per_step = per_step

    def __str__(self):
        return self.__dict__.__str__()

    def set_a_b(self, a, b):
        self.A = a
        self.B = b

    def set_vlp(self, vlp):
        self.Vlb = vlp

    def set_vnp(self, vnp):
        self.Vnb = vnp

    def calc_vwb(self):
        self.Vwb = 1 - (self.Vlb + self.Vnb)

    def set_kg(self, kg):
        self.Kg = kg

    # returns k_1 of chemical for this pplank dependent on chem kow
    def calc_k1(self, chem_kow, chem_index):
        self.k_1[chem_index] = 1/(self.A + (self.B/chem_kow))

    # returns k_2 of chemical for this pplank
    def calc_k2(self, chem_kow, chem_index, beta1, beta4, density_lip=.9, density_w=1):


        k_pw = (self.Vlb * beta1 * chem_kow) / density_lip + (self.Vnb * beta4 * chem_kow) + (self.Vwb / density_w)

        self.k_2[chem_index] = self.k_1[chem_index] / k_pw

    def init_check(self):

        checks = 0

        # check 1
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry is not None:
                count += 1

        if count == len(atts.values()):
            checks += 1

        # check 2
        if len(self.k_1) == len(self.k_2):
            checks += 1

        if checks == 2:
            return True

        return False

    def solve_steady_state(self, Cwd, i):

        return (self.k_1[i] * Cwd) / (self.k_2[i] + self.Kg)

    def solve_next_time_step(self, Cwd, chem_index, pre_step, t):
        q = self.k_1[chem_index] * Cwd
        k = (self.k_2[chem_index]) + self.Kg


        analytic = (q/k) + math.exp(-k*t) * (0-(q/k))

        top = (pre_step * ((1/self.days_per_step)-(k/2))+ q)
        bottom = (1/self.days_per_step) + (k/2)

        return analytic


class Chemical:

    def __init__(self, name, kow, cs):
        self.name = name
        self.Kow = math.pow(10, kow)
        self.Cs = cs
        self.phi = -1
        self.Cwp = -1
        self.Cwto = -1
        self.Cwdo = -1
        self.Ddoc = -1
        self.Dpoc = 0
        self.beta1 = 1
        self.beta2 = .035
        self.beta3 = 1
        self.beta4 = .35
        self.beta5 = .035
        self.Ew = self.calc_ew()
        self.Ed = self.calc_ed()

    def __str__(self):
        return self.__dict__.__str__()

    def set_cwto(self, cwto):
            self.Cwto = cwto

    def set_cwdo(self, cwdo):
            self.Cwdo = cwdo

    def set_cwp(self, cwp):

        self.Cwp = cwp

    def set_ddoc_dpoc(self, ddoc, dpoc):
        self.Ddoc = ddoc
        self.Dpoc = dpoc

    def set_beta1(self, beta1):

        self.beta1 = beta1

    def set_beta2(self, beta2):

        self.beta2 = beta2

    def set_beta3(self, beta3):

        self.beta3 = beta3

    def set_beta4(self, beta4):

        self.beta4 = beta4

    def set_beta5(self, beta5):

        self.beta5 = beta5

    def calc_ew(self):
        Ew = 1/(1.85 + (155/self.Kow))
        return Ew

    def calc_ed(self):
        Ed = math.pow((0.0000003*self.Kow + 2), -1)
        return Ed

    def calc_cs(self, region):
        
        self.Cs = region.Ocs * .35 * self.Kow * self.Cwp

    def calc_pore_water(self, Ocs):
        
        self.Cwp = self.Cs/(Ocs * .35 * self.Kow)


    def calc_phi_and_cwdo(self, region):


        adoc = region.adoc
        apoc = region.apoc
        xdoc = region.Xdoc
        xpoc = region.Xpoc

        self.phi = 1/((1+xpoc*self.Dpoc*apoc*self.Kow)+(xdoc*self.Ddoc*adoc*self.Kow))

        self.Cwdo = (self.Cwto/1000) * self.phi

        
    def init_check(self):
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry is not None:
                count += 1

        if count + 2 >= len(atts.values()):
            return True

        return False


class Region:

    def __init__(self, name, t, xdoc, xpoc, css, ocs, s, adoc=.08, apoc=.35):

        self.name = name
        self.T = t   # mean water temperture
        self.Xdoc = xdoc  # Concentration Dissolved Organic Carbon Content
        self.Xpoc = xpoc  # Concentration Particulate Organic Carbon Content
        self.Css = css  # Concentration of suspended solids in water
        self.Ocs = ocs  # Sediment organic Carbon Content
        self.S = s  # Dissolved oxygen saturation
        self.Cox = 0   # Dissolved Oxygen Concentration
        self.adoc = adoc  # DOC octanol proportionality constant
        self.apoc = apoc  # POC octanol proportionality constant

    def __str__(self):
        return str(self.__dict__)

    def set_adoc_apoc(self, apoc, adoc):
        self.adoc = adoc
        self.apoc = apoc

    def set_cox(self, cox):
            self.Cox = cox

    def calc_cox(self):
        self.Cox = (((-.24 * self.T) + 14.04) * self.S)

    def init_check(self):
        atts = self.__dict__
        count = 0
        for entry in atts.values():
            if entry is not None:
                count += 1

        if count == len(atts.values()):
            return True

        return False
