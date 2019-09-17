#################################################################
#                                                               #
#           Welcome to the ARIEL follow-up Simulator            #
#                                                               #
# This file is the main file that runs the simulator            #
# It copies an unused database file, populates with telescopes, #
# and runs the simulation                                       #
# Requires an accuracy threshold in minutes, a csv file with    #
# telescope data, and an operating mode                         #
# Generates a csv file with running total data once finished    #
#                                                               #
# Hamish Caines 07-2019                                         #
#################################################################


def parse_arguments():
    """
    Allows arguments to be parsed when this method is called. These arguments specify the setup for the simulation
    :return: threshold: accuracy threshold to be used as the cutoff in the simulation
    :return: telescope_file: name of the .csv file containing the telescopes to be used in the simulation
    :return: mode: operating mode for the simulation, controls the amount of telescope time to be used
    """
    import argparse
    parser = argparse.ArgumentParser(description='Run simulation')
    # specify arguments to be collected
    parser.add_argument('threshold', type=int, help='Accuracy threshold')  # name of file, not location
    parser.add_argument('telescopes', type=str, help='file containing telescope data for simulation')
    parser.add_argument('mode', type=str, help='Operating mode')
    args = parser.parse_args()  # collect arguments
    telescope_file = '../telescopes/' + args.telescopes  # add location of telescope file to name
    mode = args.mode
    threshold = args.threshold
    print('Target threshold: ' + str(threshold))
    print('Operating mode: ' + str(mode))
    print('Obtaining telescopes from "' + telescope_file + '"')
    return threshold, telescope_file, mode


def create_simulation_name(telescope_file, threshold, mode):
    """
    Generates folder name based on simulation settings
    :param telescope_file: name of the .csv file containing the telescopes to be used in the simulation
    :param threshold: accuracy threshold to be used as the cutoff in the simulation
    :param mode: operating mode for the simulation, controls the amount of telescope time to be used
    :return: name of folder to be used for the simulation
    """
    # stitch together name for simulation
    sim_name = telescope_file.split('/')[-1].split('.')[0] + '_' + mode + '_' + str(threshold)
    print('Simulation name: '+sim_name+', folder created')
    return sim_name+'_test'


def copy_database(sim_name):
    """
    Create folder for new simulation and copy clean database into it, renaming as new simulation name
    :param sim_name: Name of simulation
    :return: Name of database file generated
    """
    import shutil
    from os import mkdir
    database_name = sim_name+'.db'  # add database suffix
    try:
        shutil.rmtree(sim_name, ignore_errors=True)  # check for existing folder and remove
        mkdir(sim_name)
        shutil.copyfile('clean/clean.db', sim_name + '/' + database_name)  # copy clean database and rename
    except FileExistsError:
        pass

    return database_name


def populate_telescopes(database):
    """
    Obtain telescope data from csv file and store in database, and generate schedule table for each telescope
    :param database: Database object connected to database file
    """
    import database_generator
    # generate telescope data table from csv file
    telescopes = database_generator.generate_sql_table_from_csv(database.telescope_file, 'TELESCOPES', database.cursor)
    # generate schedule table for each telescope
    for telescope in telescopes:
        database.cursor.execute('DROP TABLE IF EXISTS ' + telescope[0])  # delete existing table
        database.cursor.execute('CREATE TABLE IF NOT EXISTS ' + telescope[0] + '(Target VARCHAR(25), RA DECIMAL(16,8), '
                                'Dec DECIMAL(16,8), ObsCenter DATETIME, RunStart DATETIME, RunEnd DATETIME, RunDuration'
                                ' TIME, Epoch REAL, UNIQUE(RunStart))')  # create new table
        print(telescope[0])
    print('Network has '+str(len(telescopes))+' telescopes')
    database.db.commit()


def run_simulation(database):
    """
    Run network simulation over duration of follow-up period
    :param database: Database object connected to database file
    :return counts: Array of the number of constrained targets at every forecasting period
    :return totals: Array of the total number of targets at every forecasting period
    :return dates: Array of the dates for every forecasting period
    """
    from datetime import datetime, timedelta
    interval = timedelta(days=7)  # set interval for scheduling interval
    database.sim_end = datetime(year=2030, month=6, day=12)  # end of simulation after ARIEL launch
    database.initial_expiry_calculation()
    database.sim_start = datetime.today()
    database.transit_forecast(database.sim_start, database.sim_start + interval)
    #start_date = database.find_earliest_date()  # start simulation at the date of the earliest transit in the database


    print('Running simulation from '+str(database.sim_start)+' until '+str(database.sim_end))

    #  set counter and interval for forecasting period
    time_since_forecast = timedelta(days=0)
    limit = timedelta(days=28)

    current_date = database.sim_start  # set start date

    # containers for result information
    totals = []
    counts = []
    dates = []

    while current_date < database.sim_end:  # loop while date is within simulation

        database.increment_total_night(current_date, interval)
        database.make_schedules(current_date, interval, database.mode)  # make schedules for the scheduling interval
        database.simulate_observations(current_date, interval)  # simulate the observations
        time_since_forecast += interval  # increment forecast counter
        #if time_since_forecast >= limit:  # check forecast counter

        current_date += interval  # increment date
        print('Forecasting on:', current_date)
        database.transit_forecast(current_date, current_date + limit)
        time_since_forecast = timedelta(days=0)


            # find constrained and total targets and store in arrays with date
            ####count, total = database.check_constrained(current_date)
            ###counts.append(count)
            ##totals.append(total)
            #dates.append(current_date.date())

    print('Simulation finished')


def write_count_results(counts, totals, dates, sim_name):
    """
    Write running totals of constrained and total targets, with dates to a csv file
    :param counts: Array of the number of constrained targets at every forecasting period
    :param totals: Array of the total number of targets at every forecasting period
    :param dates: Array of the dates for every forecasting period
    :param sim_name: Name of simulation
    """
    outfile = sim_name+'.csv'  # generate filename for data
    # write results to file
    with open(outfile, 'w') as f:
        f.write('# Date\tNoObserved\tNoTargets')
        for i in range(0, len(totals)):
            f.write('\n' + str(dates[i]) + '\t' + str(counts[i]) + '\t' + str(totals[i]))
    f.close()


def main():
    from os import chdir, getcwd
    import actions
    threshold, telescope_file, mode = parse_arguments()  # collect arguments
    sim_name = create_simulation_name(telescope_file, threshold, mode)  # generate simulation name
    database_name = copy_database(sim_name)  # create new database file
    chdir(sim_name)

    # create Database object connected to database file
    database = actions.Database(database_name, telescope_file, mode, threshold)
    populate_telescopes(database)

    run_simulation(database)

    #write_count_results(counts, totals, dates, sim_name)

    database.propagate_final_errors()
    count, total = database.check_constrained()
    database.store_results(count, total)


if __name__ == '__main__':
    main()
