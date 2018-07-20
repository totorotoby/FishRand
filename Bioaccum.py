# Where we run the bioaccumlation model

import Classes as obj
import FR_Input

def init_region(reg_data):


def init_phyto(phyto_data):


def init_zoop(zoo_data):


def init_chems(chem_data):



def init_fish(fish_data, diet_data, region):

    fish_log = make_fishlog(fish_data,diet_data)

    for fish in fish_log:
        # setting diet percentages
        fish.calc_diet_per(fish_log)
        # do gut percentages
        fish.calc_gut_per()
        # do Gill vent rate
        fish.calc_gv(region.Cox)
        # set feeding rate if not  already done



# returns array of fish
def make_fishlog(fish_data, diet_data):

    fish_log = []

    for i in range (len(fish_data)):
        fish = fish_data[i]
        fish_diet = diet_data.get(fish[0])
        Fish = obj.Fish(fish[0], fish[1], fish[2], fish[4], fish_diet, fish[10])
        if fish[3] != '':
            Fish.set_vnb(fish[3])
        if fish[6] != '':
            Fish.set_el(fish[6])
        if fish[7] != '':
            Fish.set_en(fish[7])
        if fish[8] != '':
            Fish.set_ew(fish[8])
        if fish[9] != '':
            Fish.set_kg(fish[9])
        Fish.set_vwb()
        fish_log.append(Fish)

    return fish_log


def main():
    reg_data, chem_data, fish_data, zoo_data, phyto_data, diet_data = FR_Input.convert_to_lists("FR_Input.xlsx")
    make_fishlog(fish_data, diet_data)
    print(fish_data)

main()