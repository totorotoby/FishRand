from main import *
from FR_Input_Output import *


data = convert_to_lists("phyto_test.xlsx", 9)
output, inputs = filter_cases(data, [9], False)
