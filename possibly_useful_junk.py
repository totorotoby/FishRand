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

    # check 2
    if len(self.k_1) == len(self.k_2):
        checks += 1

    if checks == 2:
        return True

    return False





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
        try:
            self.Gv = 1400 * ( math.pow(self.Wb, .65) / region_cox)
        except ValueError:
            print('Weight of ' + self.name + ' was sampled as a negative number. This makes no sense. Check the range of the weight distribution. \n Caught at Gv calculation.')
            exit(0)
    # sets G_d for this zooplank
    def calc_gd(self, t, css, flag,sigma=1):
        if flag == 1:
            self.Gd = self.Gv * css * sigma
        if flag == 0:
            try:
                self.Gd = .022 * math.pow(self.Wb, .85) * math.exp((.06 * t))
            except ValueError:
                print('Weight of ' + self.name + ' was sampled as a negative number. This makes no sense. Check the range of the weight distibution. \n Caught at Gd calculation.')
                exit(0)
    def calc_kg(self, T):
        if T < 15:
            try:
                self.Kg = .0005 * math.pow(self.Wb, -.2)
            except ValueError:
                print('Weight of ' + self.name + ' was sampled as a negative number. This makes no sense. Check the range of the weight distibution. \n Caught at Kg calculation.')
                exit(0)
        if T >= 15:
            try:
                self.Kg = .00251 * math.pow(self.Wb, -.2)
            except ValueError:
                print('Weight of ' + self.name + ' was sampled as a negative number. This makes no sense. Check the range of the weight distibution. \n Caught at Gv calculation.')
                exit(0)
    # returns k_1 of chemical for this zooplank
    def calc_k1(self, chem_ew):
        k1 = chem_ew * (self.Gv / self.Wb)
        return k1

    # returns k_2 of chemical for this zooplank
    def calc_k2(self, chem_kow, k_1,beta=.035, density_lip=.9 , density_w=1):

        k_bw = (self.Vlb * chem_kow)/density_lip + (self.Vnb * beta * chem_kow) + (self.Vwb/density_w)
        return k_1 / k_bw

    # returns k_d of chemical for this zooplank
    def calc_kd(self, chem_ed):

        return chem_ed * (self.Gd / self.Wb)




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