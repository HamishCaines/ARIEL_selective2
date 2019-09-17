#################################################################
# This is the main wheelhouse for the simulator                 #
# Contains the main operations needed to run the simulator      #
# Calls auxilliary files for large operations involving         #
# querying databases or mathematical functions                  #
#                                                               #
# Hamish Caines 07-2019                                         #
#################################################################


class Database:
    """
    Database object, contains database and cursor used to make changes.
    Functions within operate on the database as part of the simulation
    """
    def __init__(self, database, telescopes, mode, threshold):
        """
        Constructor that connects to, and stores information about the simulation
        :param database: name of database file for the simulation: str
        :param telescopes: location of telescope file: str
        :param mode: Operating mode to be used in the simulation: str
        :param threshold: Accuracy threshold for the simulation in minutes: int
        """
        import sqlite3
        from datetime import timedelta
        self.db = sqlite3.connect(database)
        self.cursor = self.db.cursor()
        self.names = self.read_target_names()
        self.telescope_file = telescopes  # store location
        if telescopes is not None:
            self.telescope_setup = telescopes.split('/')[-1].split('.')[0]  # store file name
        else:
            self.telescope_setup = None
        self.telescope_data = []
        self.mode = mode
        self.threshold = threshold

        self.total_night = timedelta(days=0)
        self.total_obs = timedelta(days=0)

        self.sim_start = None
        self.sim_end = None

    def update(self, table, column, row, value):
        """
        Update single value for single row in SQL table
        :param table: Table name: str
        :param column: Column name: str
        :param row: Row to update: str
        :param value: New value: int or float
        """
        self.cursor.execute('UPDATE '+table+' SET '+column+' = '+str(value)+' WHERE Name = \''+row+'\'')
        self.db.commit()

    def read_target_names(self):
        """
        Read all target names from database
        :return: List of target names
        """
        rows = self.cursor.execute('SELECT Name FROM TARGET_DATA').fetchall()
        names = []
        # extract strings from tuples
        for row in rows:
            names.append(row[0])
        return names

    def find_earliest_date(self):
        """
        Find date of earliest transit in database
        :return:
        """
        import datetime
        self.cursor.execute('SELECT Center FROM NEEDED_TRANSITS')  # find all transits
        rows = self.cursor.fetchall()
        lowest = min(rows)  # find earliest date
        self.sim_start = datetime.datetime.strptime(lowest[0], '%Y-%m-%d %H:%M:%S.%f')  # convert to datetime and store

    def run_queries(self):
        """
        Query databases to obtain known observations and missing data for each target
        Queries ETD for observations, and exoplanets.org for missing data
        """
        import query_tools
        total = len(self.names)
        count = 1
        print('Starting queries for ' + str(total) + ' targets...')
        for name in self.names:
            try:
                float(name[:-1])  # check for real target, fake star names only contain numbers
                print('Target ' + name + ' is not real (' + str(count) + '/' + str(total) + ')')
            except ValueError:  # catch error for real targets
                print('Querying ' + name + ' (' + str(count) + '/' + str(total) + ')')
                target = query_tools.Target(name)  # initialise Target object
                target.ETD_query()  # obtain observation data
                # check for missing data and query exoplanets.org to fill gaps
                if list(target.__dict__.values()).__contains__(None):
                    target.EXO_query()
                    target.EXO_used = True
                target.write_query_data(self.cursor, self.db)  # write results
            count += 1

        print('Initial data obtained, now fitting periods...')
        self.initial_period_fit()  # perform initial period fit for initial observation available
        # self.initial_prop_to_ariel()  # calculate uncertainties at ARIEL launch based on current observation data

    def initial_period_fit(self):
        """
        Perform initial period fit for each target based on initial observation data
        :return:
        """
        import data_tools
        total = len(self.names)
        count = 1
        for target in self.names:
            obs = data_tools.read_obs_data(self.cursor, target)  # load target data
            try:
                print('Fitting ' + target + ' ('+str(count)+'/'+str(total)+')')
                # perform fit
                fit_period, fit_period_err, latest_tmid, latest_tmid_err, latest_epoch = data_tools.period_fit(obs)
                # store values in table
                self.update('TARGET_DATA', 'FitPeriod', target, fit_period)
                self.update('TARGET_DATA', 'FitPeriodErr', target, fit_period_err)
                self.update('TARGET_DATA', 'CurrentPeriod', target, fit_period)
                self.update('TARGET_DATA', 'CurrentPeriodErr', target, fit_period_err)
                self.update('TARGET_DATA', 'TruePeriod', target, fit_period)
                self.update('TARGET_DATA', 'TruePeriodErr', target, fit_period_err)
                self.update('TARGET_DATA', 'TrueEpoch', target, latest_epoch)
                self.update('TARGET_DATA', 'TrueLastObs', target, latest_tmid)
                self.update('TARGET_DATA', 'TrueLastObsErr', target, latest_tmid_err)

                # check for missing error data on last observation and store latest value with error available
                if self.cursor.execute('SELECT LastObsErr FROM TARGET_DATA WHERE Name = \''+target+'\'').fetchall()[0]\
                        is None:
                    self.update('TARGET_DATA', 'LastObs', target, latest_tmid)
                    self.update('TARGET_DATA', 'LastObsErr', target, latest_tmid_err)
                    self.update('TARGET_DATA', 'LastEpoch', target, latest_epoch)

                # check for new latest observation
                elif latest_tmid > self.cursor.execute('SELECT LastObs FROM TARGET_DATA WHERE '
                                                       'Name = \''+target+'\'').fetchall()[0][0]:
                    self.update('TARGET_DATA', 'LastObs', target, latest_tmid)
                    self.update('TARGET_DATA', 'LastObsErr', target, latest_tmid_err)
                    self.update('TARGET_DATA', 'LastEpoch', target, latest_epoch)

            except Warning:  # thrown by period_fit
                print('Fit for ' + target + ' failed, has '+str(len(obs))+' observations')
                self.cursor.execute(
                    'UPDATE TARGET_DATA SET CurrentPeriod = PeriodStart, CurrentPeriodErr = PeriodStartErr WHERE '
                    'FitPeriod IS NULL')
                self.cursor.execute(
                    'UPDATE TARGET_DATA SET TruePeriod = PeriodStart, TruePeriodErr = PeriodStartErr WHERE FitPeriod IS'
                    ' NULL')
                self.cursor.execute('UPDATE TARGET_DATA SET TrueEpoch = LastEpoch WHERE TrueEpoch IS NULL')

            finally:
                count += 1
                self.db.commit()

    def initial_prop_to_ariel(self):
        """
        Calculate timing uncertainty at ARIEL launch for each target based on current observation data available
        """
        import data_tools

        # obtain required data for all targets
        self.cursor.execute(
            'SELECT Name, CurrentPeriod, CurrentPeriodErr, LastObs, LastObsErr FROM TARGET_DATA '
            'WHERE CurrentPeriodErr NOT NULL AND LastObsErr NOT NULL')
        rows = self.cursor.fetchall()
        # for each target
        for row in rows:
            name = row[0]
            try:
                err_tot, percent, loss = data_tools.prop_to_date(row, self.sim_end)  # calculate propagated error
                # store values in table
                self.update('TARGET_DATA', 'ErrAtAriel', name, err_tot)
                self.update('TARGET_DATA', 'PercentLoss', name, percent)
                self.update('TARGET_DATA', 'LossAtAriel', name, loss)
                self.update('TARGET_DATA', 'ErrAtArielStart', name, err_tot)
                self.update('TARGET_DATA', 'PercentLossStart', name, percent)
                self.update('TARGET_DATA', 'LossAtArielStart', name, loss)
            except Warning:
                pass

    def initial_expiry_calculation(self):
        import data_tools
        #self.cursor.execute('ALTER TABLE TARGET_DATA ADD COLUMN IF NOT EXISTS ExpiryJD DECIMAL(16,8)')
        self.db.commit()
        rows = self.cursor.execute(
            'SELECT Name, CurrentPeriod, CurrentPeriodErr, LastObs, LastObsErr, ExpiryJD FROM TARGET_DATA WHERE CurrentPeriodErr N'
            'OT NULL AND LastObsErr NOT NULL AND Depth > 10.0').fetchall()
        for row in rows:
            name = row[0]
            expiry_jd = data_tools.prop_to_threshold(row, self.threshold, self.sim_end)

            self.update('TARGET_DATA', 'ExpiryJD', name, expiry_jd)
        self.db.commit()

    def transit_forecast(self, start, end):
        """
        Forecast exoplanet transits visible from somewhere on earth within the window specified by arguments
        :param start: Start of forecast window: datetime
        :param end: End of forecast window: datetime
        """
        import observation_tools
        import julian

        #  make table for transit data if not exists
        self.cursor.execute(
            'CREATE TABLE IF NOT EXISTS NEEDED_TRANSITS(Center DATETIME, Name VARCHAR(20), Ingress DATETIME, Egress DATET'
            'IME, Duration TIME, RA DECIMAL(9,7), Dec DECIMAL(9,7), PercentLoss REAL, Epoch REAL, ErrAtAriel DECIMAL(16'
            ',8), ExpiryJD DECIMAL(16,8))')

        start_jd = julian.to_jd(start, fmt='jd') - 2400000  # convert start to JD format in table
        # select targets with deep transits and discovery dates in the past, AND require further observation
        #print('SELECT Name FROM TARGET_DATA WHERE Depth > 10.0 AND ExpiryJD < ' + str(start_jd))

        rows = self.cursor.execute('SELECT Name, ExpiryJD FROM TARGET_DATA WHERE Depth > 10.0 AND ExpiryJD < ' + str(start_jd)).fetchall()
        print('needed', rows)

        # obtain deep names and reset names in database
        deep_names = []
        for row in rows:
            deep_names.append(row[0])
        self.names = deep_names

        # for each deep target
        for name in self.names:
            data = observation_tools.read_data(self.cursor, name)
            transits = observation_tools.transit_forecast(data, name, start, end)  # forecast transits for target

            for transit in transits:
                # add new transit to the table
                string = 'INSERT INTO NEEDED_TRANSITS (Center, Name, Ingress, Egress, Duration, RA, Dec, PercentLoss, ' \
                         'Epoch, ErrAtAriel) VALUES ( \'' + str(transit.center) + '\', \'' + transit.name + '\', \'' \
                         + str(transit.ingress) + '\', \'' + str(transit.egress) + '\', \'' + str(
                            transit.duration) + '\', ' + str(transit.ra) + ', ' + str(transit.dec) + ', '
                # check for missing error values and modify string accordingly
                if transit.error is None:
                    string += 'NULL, ' + str(transit.epoch) + ', NULL)'
                else:
                    string += str(transit.loss) + ', ' + str(transit.epoch) + ', ' + str(transit.error) + ')'
                self.cursor.execute(string)

        self.db.commit()

    def increment_total_night(self, start, interval):
        """
        Keep a running total of the total available observing hours through out simulation by calculating sunset and
        rise for each day in specified window
        :param start: Start of window: datetime
        :param interval: Length of window: datetime
        """
        from datetime import timedelta
        import mini_staralt

        end = start + interval
        day = timedelta(days=1)

        while start < end:
            for telescope in self.telescope_data:  # calculate sunset/rise at each site
                sunset, sunrise = mini_staralt.sun_set_rise(start, lon=telescope.lon, lat=telescope.lat,
                                                            sundown=-12)
                duration = sunrise - sunset
                self.total_night += duration  # add duration to counter

            start += day  # increment day

    def load_telescope_data(self):
        """
        Create Telescope objects for telescopes in database and store in Database object
        """
        import observation_tools
        self.cursor.execute('SELECT * FROM TELESCOPES')  # query database
        rows = self.cursor.fetchall()
        telescopes = []
        # for each telescope
        for row in rows:
            new_telescope = observation_tools.Telescope()
            new_telescope.gen_from_database(row)
            telescopes.append(new_telescope)

        self.telescope_data = telescopes

    def obtain_visible_transits(self, start, interval):
        import observation_tools
        end = start + interval

        candidates = self.cursor.execute(
            'SELECT * FROM NEEDED_TRANSITS WHERE Center BETWEEN "' + str(start.date()) + ' 00:00:00" AND "' + str(
                end.date()) + ' 00:00:00" ORDER BY Center ASC').fetchall()
        self.load_telescope_data()
        transits = []
        for candidate in candidates:
            new_transit = observation_tools.Transit()
            new_transit.gen_from_database(candidate)
            if new_transit.check_visibility_telescopes(self.telescope_data):
                transits.append(new_transit)

        return transits

    def make_schedules(self, start, interval, mode):
        import numpy as np

        transits = self.obtain_visible_transits(start, interval)
        end = start + interval

        for telescope in self.telescope_data:
            matching_transits = []
            for transit in transits:
                if telescope.name in transit.telescope:
                    matching_transits.append(transit)

            if mode == 'unlimited':
                limit = np.inf
            elif 'perweek' in mode:
                limit = int(mode[0])
            else:
                print('NO VALID MODE SPECIFIED')
                raise IOError

            self.schedule(matching_transits, limit, telescope.name, end)

    def schedule(self, transits, limit, telescope, end):
        """
        Make schedule for a given telescope given a set of transits and a limit on the number of observations
        available in the interval
        :param transits: List of Transit objects for transits in the given window visible from the telescope
        :param limit: Maximum number of observations allowed in the interval
        :param telescope: Name of telescope being scheduled
        """
        from datetime import datetime, timedelta
        from sqlite3 import IntegrityError

        new_transits = 0  # set limit counter
        while new_transits < limit:  # loop while limit not reached
            # iterate through transits
            for transit in transits:
                # check already scheduled observations
                scheduled = self.cursor.execute('SELECT RunStart, RunEnd FROM "' + telescope + '" WHERE RunEnd < "'+str(end)+'"').fetchall()
                # set start and finish times including continuum observation time
                continuum = timedelta(minutes=45)
                new_start = transit.ingress - continuum
                new_end = transit.egress + continuum
                duration = new_end - new_start  # find duration
                # check for empty schedule and schedule transit if empty
                if len(scheduled) == 0:
                    try:
                        self.cursor.execute(
                            'INSERT INTO ' + telescope + ' VALUES ("' + transit.name + '", ' + str(transit.ra) +
                            ', ' + str(transit.dec) + ', "' + str(transit.center) + '", "' + str(
                                new_start) + '", "' + str(new_end) + '", "' +
                            str(duration) + '", ' + str(transit.epoch) + ')')

                        new_transits += 1  # increment number of scheduled observations
                        self.total_obs += transit.duration  # add to total scheduled time
                    except IntegrityError:
                        pass
                # if schedule is not empty, need to check for space in schedule
                else:
                    space = True
                    for single in scheduled:  # check for each scheduled observations
                        # obtain start and end times
                        try:
                            start_time = datetime.strptime(single[0], '%Y-%m-%d %H:%M:%S.%f')
                            end_time = datetime.strptime(single[1], '%Y-%m-%d %H:%M:%S.%f')
                        except ValueError:
                            start_time = datetime.strptime(single[0], '%Y-%m-%d %H:%M:%S')
                            end_time = datetime.strptime(single[1], '%Y-%m-%d %H:%M:%S')

                        if start_time < new_start < end_time:  # new observation starts before current one ends
                            space = False
                        if start_time < new_end < end_time:  # new observation ends after current one starts
                            space = False
                        if new_start < start_time and end_time < new_end:  # new observation surrounds current one
                            space = False

                    if space:  # schedule if space is available
                        try:
                            self.cursor.execute(
                                'INSERT INTO ' + telescope + ' VALUES ("' + transit.name + '", ' + str(
                                    transit.ra) +
                                ', ' + str(transit.dec) + ', "' + str(transit.center) + '", "' + str(
                                    new_start) + '", "' + str(new_end) + '", "' +
                                str(duration) + '", ' + str(transit.epoch) + ')')

                            new_transits += 1  # increment number of scheduled observations
                            self.total_obs += transit.duration  # add to total scheduled time
                        except IntegrityError:
                            pass


            # break when limit reached
            break

    def simulate_observations(self, start, interval):
        import observation_tools
        end = start + interval
        counter = 0
        for telescope in self.telescope_data:
            observations = self.cursor.execute('SELECT * FROM "' + telescope.name + '" WHERE ObsCenter BETWEEN "' + str(
                start.date()) + ' 00:00:00" AND "' + str(end.date()) + ' 00:00:00"').fetchall()
            for observation in observations:
                success = observation_tools.flip_unfair_coin()
                if success:
                    counter += 1
                    data = self.cursor.execute(
                        'SELECT LastObs, LastEpoch, CurrentPeriod, TrueLastObs, TruePeriod,'
                        ' TrueEpoch FROM TARGET_DATA WHERE Name = "' + observation[0] + '"').fetchall()[0]
                    new_epoch = observation[-1]
                    last_tmid, last_epoch, period, true_last_tmid, true_period, true_epoch = \
                        data[0], data[1], data[2], data[3], data[4], data[5]

                    new_tmid, new_tmid_err = observation_tools.generate_results(last_tmid, last_epoch, new_epoch,
                                                                                period)
                    exp = self.cursor.execute('SELECT ExpiryJD FROM TARGET_DATA WHERE Name = "'+observation[0]+'"').fetchall()
                    if exp[0][0] > new_tmid:
                        print(exp, new_tmid)
                        print()
                    if true_last_tmid is not None:
                        true_center = observation_tools.find_true_t0(true_last_tmid, true_period, true_epoch, new_epoch)
                    else:
                        true_center = 'NULL'
                    # add new observation to the database
                    self.add_new_observation(observation[0], new_epoch, new_tmid, new_tmid_err,
                                             telescope.name, true_center)
                    self.db.commit()

        print('observed', counter, 'from', start, 'to', end)

    def add_new_observation(self, name, epoch, tmid, tmid_err, telescope, true_center):
        """
        Add a new observation to the database and trigger recalculation period and error at ARIEL launch
        :param name: Name of target: str
        :param epoch: Epoch of observation: int
        :param tmid: Observed transit center: float
        :param tmid_err: Error in observed transit center: float
        :param telescope: Name of telescope used: str
        :param true_center: Expected ephemeris based on expected period: float
        :return:
        """
        import observation_tools
        from sqlite3 import IntegrityError

        self.update('TARGET_DATA', 'LastObs', name, tmid)
        self.update('TARGET_DATA', 'LastObsErr', name, tmid_err)
        self.update('TARGET_DATA', 'LastEpoch', name, epoch)

        # find new observation ID
        ids_raw = self.cursor.execute('SELECT ObID FROM "' + name + '"').fetchall()
        new_id = observation_tools.find_highest_id(ids_raw)

        # insert into database
        try:
            self.cursor.execute('INSERT INTO "' + name + '" (ObID, Epoch, ObsCenter, ObsCenterErr, TrueCenter, Source) '
                                                         'VALUES (' + str(new_id) + ', ' + str(epoch) + ', '
                                + str(tmid) + ', ' + str(tmid_err) + ', ' + str(true_center) + ', "' + telescope + '")')
            self.cursor.execute('UPDATE TARGET_DATA SET NoOfObs = NoOfObs + 1 WHERE Name = "' + name + '"')
            self.db.commit()
        except IntegrityError:
            pass

        self.recalculate(name)  # trigger recalculation based on new data

    def recalculate(self, name):
        """
        Recalculate period and error at ARIEL launch based on data currently in database
        :param name: Target name to be recalculated: str
        """
        import data_tools

        ######################
        # PERIOD FIT SECTION #
        ######################

        # obtain current observations
        observations = data_tools.read_obs_data(self.cursor, name)

        # fit new period and update table values
        try:
            new_period, new_period_err, tmid_max, tmid_max_err, last_epoch = data_tools.period_fit(observations)
            self.update('TARGET_DATA', 'CurrentPeriod', name, new_period)
            self.update('TARGET_DATA', 'CurrentPeriodErr', name, new_period_err)
            self.update('TARGET_DATA', 'LastObs', name, tmid_max)
            self.update('TARGET_DATA', 'LastObsErr', name, tmid_max_err)
            self.update('TARGET_DATA', 'LastEpoch', name, last_epoch)
            self.db.commit()

            if observations[-1].tmid_err < self.threshold / 24 / 60:
                row = self.cursor.execute(
                    'SELECT Name, CurrentPeriod, CurrentPeriodErr, LastObs, LastObsErr, ExpiryJD FROM TARGET_DATA WHERE CurrentPeriodErr N'
                    'OT NULL AND LastObsErr NOT NULL AND Depth > 10.0 AND Name = "' + name + '"').fetchall()[0]
                expiry_jd = data_tools.prop_to_threshold(row, self.threshold, self.sim_end)
                if name == 'WASP-74b':
                    print('hello')
                self.update('TARGET_DATA', 'ExpiryJD', name, expiry_jd)
            #self.update('TARGET_DATA', 'TrueLastObs', name, true_center)

        except Warning:
            # insufficient observations for period fit, 3 required
            if observations is None:
                print('Fit for '+name+'failed, has no observations')
            else:
                if len(observations) != 1:
                    print('Fit for '+name+' failed, has '+str(len(observations))+' observations')
                else:
                    print('Fit for ' + name + ' failed, has ' + str(len(observations)) + ' observation')
        except TypeError:
            pass
        self.db.commit()

    def store_results(self, count, total):
        """
        Store running total data in results table
        :param count: Number of constrained targets: int
        :param total: Total targets in simulation: int
        """
        import sqlite3
        # create results table
        self.cursor.execute('CREATE TABLE IF NOT EXISTS RESULTS (Network VARCHAR(25), Mode VARCHAR(25), Accuracy REAL, '
                            'Constrained REAL, Total REAL, PercentTargets REAL, TotalNight DATETIME, TotalObsTime DATETIME, PercentUsage REAL)')
        percent_constrained = count/total*100

        # convert running night and observation time totals into seconds
        total_night_s = self.total_night.total_seconds()
        total_obs_s = self.total_obs.total_seconds()
        percent_usage = total_obs_s/total_night_s*100  # calculate percentage

        # store values
        self.cursor.execute('INSERT INTO RESULTS VALUES ("'+self.telescope_setup+'", "'+self.mode+'", '
                            ''+str(self.threshold)+', '+str(count)+', '+str(total)+', '+str(percent_constrained)+', "'
                            + str(self.total_night)+'", "'+str(self.total_obs)+'", '+str(percent_usage)+')')
        self.db.commit()

        results_db = sqlite3.connect('../results.db')
        results_cursor = results_db.cursor()
        results_cursor.execute('CREATE TABLE IF NOT EXISTS RESULTS (Network VARCHAR(25), Mode VARCHAR(25), Accuracy REAL, '
                            'Constrained REAL, Total REAL, PercentTargets REAL, TotalNight TIME, TotalObsTime TIME, PercentUsage REAL)')
        results_cursor.execute('INSERT INTO RESULTS VALUES ("'+self.telescope_setup+'", "'+self.mode+'", '
                            ''+str(self.threshold)+', '+str(count)+', '+str(total)+', '+str(percent_constrained)+', "'
                            + str(self.total_night)+'", "'+str(self.total_obs)+'", '+str(percent_usage)+')')
        results_db.commit()

        with open('../results.csv', 'a+') as f:
            f.write(self.telescope_setup+'", "'+self.mode+'", '
                            ''+str(self.threshold)+', '+str(count)+', '+str(total)+', '+str(percent_constrained)+', "'
                            + str(self.total_night)+'", "'+str(self.total_obs)+'", '+str(percent_usage)+'\n')
            f.close()

    def check_constrained(self):
        names = self.cursor.execute('SELECT Name FROM TARGET_DATA WHERE Depth > 10.0').fetchall()
        total = len(names)
        constrained = 0
        for name in names:
            err_at_ariel = self.cursor.execute('SELECT ErrAtAriel FROM TARGET_DATA WHERE Name = "'+name[0]+'"').fetchall()[0]
            if err_at_ariel[0] < self.threshold/24/60:
                constrained += 1

        return constrained, total

    def propagate_final_errors(self):
        import julian
        names = self.cursor.execute('SELECT Name FROM TARGET_DATA WHERE Depth > 10.0').fetchall()
        for name in names:
            data = self.cursor.execute('SELECT CurrentPeriod, CurrentPeriodErr, LastObs, LastObsErr FROM TARGET_DATA WHERE Name = "'+name[0]+'" AND Depth > 10.0').fetchall()
            if data[0][2] > julian.to_jd(self.sim_end, fmt='jd') - 2400000:
                self.update('TARGET_DATA', 'ErrAtAriel', name[0], data[0][3])
            else:
                err_at_ariel = self.prop_to_date(data[0], self.sim_end)
                self.update('TARGET_DATA', 'ErrAtAriel', name[0], err_at_ariel)

            self.db.commit()

    def prop_to_date(self, row, sim_end):
        import julian
        import numpy as np

        end_jd = julian.to_jd(sim_end, fmt='jd') - 2400000

        period = row[0]
        period_err = row[1]
        tmid = row[2]
        tmid_err = row[3]

        current = tmid
        epochs = 0

        if tmid_err is not None and period_err is not None:
            if current < end_jd:
                while current < end_jd:
                    epochs += 1
                    current += period

                err_at_ariel = np.sqrt(tmid_err*tmid_err + epochs*epochs*period_err*period_err)

            else:
                err_at_ariel = tmid_err
        else:
            err_at_ariel = 1000

        return err_at_ariel
