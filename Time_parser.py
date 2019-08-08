from datetime import datetime

def num_steps(start_string, end_string, step):

    start = datetime.strptime(start_string, '%m,%d,%Y').date()
    end = datetime.strptime(end_string, '%m,%d,%Y').date()

    if step == 'Weeks':
        difference = (end-start).days//7
        time_per_step = 7
    elif step == 'Months':
        difference = (end-start).days//30
        time_per_step = 30
    elif step == 'Quarters':
        difference = (end - start).days // 91
        time_per_step = 91
    elif step == 'Days':
        difference = (end - start).days - 1
        time_per_step = 1
    else:
        print('Check that your choice of timesteps is avaliable.')
        exit(0)
    return difference, time_per_step
