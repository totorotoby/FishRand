from datetime import datetime

#TODO make output not just in time steps but in date on excel tabs
def num_steps(start_string, end_string, step):

    start = datetime.strptime(start_string, '%m,%d,%Y').date()
    end = datetime.strptime(end_string, '%m,%d,%Y').date()

    if step == 'Weeks':
        difference = (end-start).days//7
        time_per_step = 7
    if step == 'Months':
        difference = (end-start).days//30
        time_per_step = 30
    if step == 'Quarters':
        difference = (end - start).days // 91
        time_per_step = 91
    if step == 'Days':
        difference = (end - start).days - 1
        time_per_step = 1

    return difference, time_per_step
