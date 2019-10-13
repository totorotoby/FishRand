# ParameterCheck.py
# Script to check that parameters match between python inputs and Aquaweb inputs

import xlrd
import sys

def openWorkbooks(pFile, aFile):

    python = xlrd.open_workbook(pFile)
    aquaWeb = xlrd.open_workbook(aFile)

    return python, aquaWeb


def compRegions(pRegion, pT, aRegion):

    print('REGIONAL')
    names = ['Temp', 'Con Dis Carbon Content', 'Con Part Carbon Content',
             'Con Sup solds', 'Frac Carbon Content in Sed', 'Dis Ox Sat',
             'Dis Ox Con']
    aquaValues = [aRegion.col(2)[i].value for i in range(11, 19)]
    aquaValues.pop(5)
    pythonValues = [pRegion.col(1)[i].value for i in range(2, 8)]
    pTemp = [pT.col(1)[1].value]
    pythonValues = pTemp + pythonValues
    for i in range(len(aquaValues)):
        if pythonValues[i] != aquaValues[i]:
            print(names[i] + ' is miss-matched')
            print( 'with ' + 'p: ' + str(pythonValues[i]) + ' and a: ' + str(aquaValues[i]))
            if pythonValues[i] == '':
                print('and  is not entered in python sheet')
            if aquaValues[i] == '':
                print('and is not entered in aquaweb sheet')
    
    print('--------------------------------------------')


def get_data(entry_col, instance_len, new_list):

    new_entry = []
    for i in range(len(entry_col)):
        if i % (instance_len + 1) == 0:
            entry = entry_col[i].value
            if len(new_entry) != 0:
                new_list.append(new_entry)
            new_entry = []
            if entry == 'END':
                break
        else:
            if i % (instance_len + 1) == 1:
                if entry_col[i].value == '':
                    break
                else:
                    new_entry.append(entry_col[i].value)
            else:
                new_entry.append(entry_col[i].value)
    if len(new_entry) != 0:
        new_list.append(new_entry)



def get_con_data(con_sheet, con_data, num_reg, num_chem, timesteps):

    for i in range (1, timesteps+1):
        t_cons = []
        column = con_sheet.col(1)
        for j in range (num_reg):
            r_cons = []
            for k in range (0, num_chem):
                c_cons = []
                for p in range (1,5):
                        
                    # This is a mess so...
                    # num_chem * 4 * j counts all red spaces in all of the  prior regional sections
                    # 1 + j counts the number of yellow regional labels
                    # k*4 counts the number of red spaces from this region to the partiular chemical
                    # 2*j + 1 + k counts the number of green labels for chemicals
                    # p counts the remaining red spaces down to the type of concentration to append
                        
                        
                    c_cons.append(column[(1+j) + ((j * num_chem) + (1+k)) + ((j * num_chem * 4) + (k * 4) + p)].value)
                r_cons.append(c_cons)
            t_cons.append(r_cons)
        con_data.append(t_cons)


def compChems(pythonChemProp, pythonChemCon, aquaChems):

    
    print('CHEMICALS')
    names = ['Name', 'Log Octanol-Water Partition Coefficent', 'Disequilbrium factor of Disolved Organic Carbon', 'Disequilbrium factor of Particulate Organic Carbon', 'Beta Lipid Omn', 'Beta Non-lipid OM', 'Beta Lipid digesta', 'Beta Sediment OC', 'Beta Non-lipid digesta', 'Concentration in Sediment', 'Total Concentration in Water', 'Dissolved Concentration in Water', 'Concentration in Pore Water']
    pythonChemPropL = []
    pythonChemConL = []
    get_data(pythonChemProp.col(1), 9, pythonChemPropL)
    numChem = len(pythonChemPropL)
    get_con_data(pythonChemCon, pythonChemConL, 1, numChem, 1)
    pythonChemConL = pythonChemConL[0][0]
    allPython = []
    for i in range(numChem):
        allPython.append(pythonChemPropL[i] + pythonChemConL[i])
        
    #getting aquaWeb inputs
    awChems = []
    for i in range(9, numChem + 9):
        awChems.append([aquaChems.row(i)[j].value for j in range(3,6)])

    for i in range(numChem):
        for j in range(len(allPython[0])):
            if j == 1:
                if allPython[i][j] != awChems[i][0]:
                    print('In ' + allPython[i][0]+ ' ' + names[j] + ' is miss-matched')
                    if type(allPython[i][j]) == str:
                        print('with p: ' + allPython[i][j] + ' and a: ' + str(awChems[i][0]))
                print()
            
            if j in [2, 3, 4, 6] and allPython[i][j] != 1 and allPython[i][j] != '':
                 print('In ' + allPython[i][0]+ ' ' + names[j] + ' is miss-matched')
                 print()
            if j in [5,8] and allPython[i][j] != .035 and allPython[i][j] != '':
                 print('In ' + allPython[i][0]+ ' ' + names[j] + ' is miss-matched')
                 print()

            if j in [2,3] and allPython[i][j] == '':
                print('In ' +allPython[i][0]+ ' ' + names[j] + ' is not entered in Python')
                print()

            if j == 9 and allPython[i][j] != awChems[i][2]:
                print('In ' + allPython[i][0] + ' ' + names[j] + ' is miss-matched')
                if type(allPython[i][j]) == str and allPython[i][j] != '':
                    print('with p: ' + allPython[i][j] + ' and a: ' + str(awChems[i][2]))
                print()
                        
            if j == 10 and allPython[i][j] != awChems[i][1]/1000:
                print('In ' + allPython[i][0] + ' ' + names[j] + ' is miss-matched')
                if type(allPython[i][j])  and allPython[i][j] != '':
                     print('with p: ' + allPython[i][j] + ' and a: ' + str(awChems[i][1]))
                print()


def compLowerOrgs(pythonOrgs, aquaOrgs):
    print('---------------------------------------------------')
    print('PHYTOS')
    pPhyto = []
    pNames = ['Growth Rate', 'Lipid Content', 'Non-Lipid Content']
    get_data(pythonOrgs.col(7), 4, pPhyto)
    pPhyto = pPhyto[0][1:]
    aquaCol = aquaOrgs.col(2)
    aPhyto = [aquaCol[27].value] +  [aquaCol[i].value for i in range(23, 25)]

    for i in range(len(aPhyto)):
        if round(aPhyto[i],2) != pPhyto[i]:
            if pPhyto[i] == '':
                pass
            else:
                print(pNames[i] + ' is miss-matched')
                if type(pPhyto[i]) == str:
                    print('with p: ' + pPhyto[i] + ' and a: ' + str(aPhyto[i]))
                print()

    print('---------------------------------------------------')
    print('ZOOPS')
    zNames = ['Weight', 'Lipid Content', 'Non-Lipid Content', 'Fraction Pore Water Ventilated', 'Dietary absorption efficency of lipid Organic Matter', 'Dietary absorption efficency of nonlipid Organic Matter', 'Dietary absorption efficency of water', 'Feeding Rate', 'Growth Rate', 'Filter Feeder Flag']
    pZoop = []
    get_data(pythonOrgs.col(5), 11, pZoop)
    pZoop = pZoop[0][1:8] + [pZoop[0][len(pZoop[0])-1]]
    aZoop = [aquaCol[i].value for i in range(33,36)] + [aquaCol[38].value] + [aquaCol[i].value for i in range(45,48)] + [aquaOrgs.col(3)[31].value]

    for i in range(len(aZoop)):
        if aZoop[i] != pZoop[i]:
            if i in range(2,9) and pZoop[i] == '':
                pass
            else:
                print(zNames[i] + ' is miss-matched')
                if type(pZoop[i]) == str:
                    print('with p: ' + pZoop[i] + ' and a: ' + str(aZoop[i]))
                print()


    
def compFish(pSheetOrg, pSheetDiet, aquaSheet):

    print('---------------------------------------------------')
    print('INVERTEBRATES')

    
    #getting python data
    pInvert = []
    get_data(pSheetOrg.col(3), 11, pInvert)
    for i in range(len(pInvert)):
        pInvert[i] = pInvert[i][:8] + [pInvert[i][len(pInvert[0])-1]]
    #getting aqua data
    aquaCol = aquaSheet.col(2) 
    numInvert = len(pInvert)
    aInvert = getAquaFish(aquaCol, aquaSheet.col(3), 56, numInvert)
    fishPrint(aInvert, pInvert, numInvert, 0)


    print('---------------------------------------------------')
    print('FISH')
    print()
    #getting python data
    pFish = []
    get_data(pSheetOrg.col(1), 11, pFish)
    for i in range(len(pFish)):
        pFish[i] = pFish[i][:8] + [pFish[i][len(pFish[0])-1]]
    #getting aqua data
    aquaCol = aquaSheet.col(2) 
    numFish = len(pFish)
    aFish = getAquaFish(aquaCol, aquaSheet.col(3), 176, numFish)
    fishPrint(aFish, pFish, numFish, 1)
                        
def fishPrint(aInvert, pInvert, numInvert, flag):

    zNames = ['Weight', 'Lipid Content', 'Non-Lipid Content', 'Fraction Pore Water Ventilated', 'Dietary absorption efficency of lipid Organic Matter', 'Dietary absorption efficency of nonlipid Organic Matter', 'Dietary absorption efficency of water', 'Feeding Rate', 'Growth Rate', 'Filter Feeder Flag']
    if flag == 0:
        innerloop = len(aInvert[0])
    else:
        innerloop = len(aInvert[0])-1
    
    for i in range(numInvert):
        for j in range(innerloop):
            if round(aInvert[i][j], 6) != pInvert[i][j+1]: 
                if j in range(2,9) and pInvert[i][j+1] == '':
                    pass
                else:
                    print('In '+ pInvert[i][0] +' ' + zNames[j] + ' is miss-matched')
                    if type(pInvert[i][j+1]) == str:
                        print('with p: ' + pInvert[i][j+1] + ' and a: ' + str(aInvert[i][j]))
                    print()

    
    
    

def getAquaFish(col1, col2, startNum, numOrg):


    data = []
    for i in range(numOrg):
        invert = []
        bIndex = startNum+(i*23) + i
        for j in range(1, 4):
            invert.append(col1[bIndex + j].value)
        invert.append(col1[bIndex + 6].value)
        for j in range(13, 16):
            invert.append(col1[bIndex + j].value)
        invert.append(col2[bIndex-1].value)
        data.append(invert)
        
        
    return data



                
def main():

    pFile = sys.argv[1]
    aFile = sys.argv[2]

    python, aquaWeb = openWorkbooks(pFile, aFile)
    print('\n\n')
    compRegions(python.sheet_by_index(1), python.sheet_by_index(2), aquaWeb.sheet_by_index(3))

    compChems(python.sheet_by_index(3), python.sheet_by_index(4), aquaWeb.sheet_by_index(5))

    compLowerOrgs(python.sheet_by_index(5), aquaWeb.sheet_by_index(4))
    compFish(python.sheet_by_index(5), python.sheet_by_index(6), aquaWeb.sheet_by_index(4))
    
main()
