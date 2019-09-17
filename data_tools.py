#################################################################
# Code to fit periods and calculate propagated uncertainties    #
# based on sets of observation data, either from ETD or the     #
# simulator                                                     #
# Called by methods in query_tools.py and actions.py            #
#                                                               #
# Hamish Caines 07-2019                                         #
#################################################################


def read_obs_data(cursor, name):
    """
    Read all observations for a given target and store as a list of Observation objects
    :param cursor: Cursor connected to database: cursor
    :param name: Target to be queried: str
    :return: List of Observation objects
    """
    # obtain observation data
    cursor.execute('SELECT Epoch, ObsCenter, ObsCenterErr FROM \''+name+'\'')
    rows = cursor.fetchall()
    observations = []
    # loop though observations
    for row in rows:
        # create Observation for each complete set of data and add to list
        if None not in row:
            new_ob = Observation(row)
            observations.append(new_ob)
    return observations


def period_fit(observations):
    """
    Performs a fit period for a given set of observations of a target
    :param observations: List of Observation objects
    :return fit_period: Resulting period from fit: float
    :return fit_period_err: Error in resulting period from fit: float
    :return tmid_max: Latest tmid in set of observations: float
    :return tmid_max_err: Error in latest tmid in set of observations: float
    :return max(epochs): Transit epoch for the latest tmid in sample
    """
    import numpy as np
    # initialise containers for values to be fitted
    epochs = []
    tmids = []
    weights = []
    # initialise tmid and err containers
    tmid_max = 0
    tmid_max_err = 0
    true_tmid_max = 0
    # loop through observations
    for ob in observations:
        # extract data
        epochs.append(ob.epoch)
        tmids.append(ob.tmid)
        weights.append(ob.weight)

        # check and set latest tmid w/error
        if ob.tmid > tmid_max:
            tmid_max = ob.tmid
            tmid_max_err = ob.tmid_err

    # attempt polynomial fit
    try:
        # check for sufficient values
        if len(epochs) >= 3:
            # run fit and extract results
            poly_both = np.polyfit(epochs, tmids, 1, cov=True, w=weights)

            poly, cov = poly_both[0], poly_both[1]
            fit_period = poly[0]
            fit_period_err = np.sqrt(cov[0][0])  # calculate error

            return fit_period, fit_period_err, tmid_max, tmid_max_err, max(epochs) #, true_tmid_max
        else:
            # if insufficient observations, raise warning
            raise Warning

    except ValueError:
        pass
    except np.linalg.LinAlgError:
        pass
    except TypeError:
        pass


def check_better_period(start, fit, fit_err):
    dif = abs(start - fit)
    if fit_err < dif:
        return True
    else:
        return False


class Observation:
    """
    Observation object, contains data for a single observation for a target
    Used in period_fit and prop_forwards
    """
    def __init__(self, row):
        """
        Constructor, extracts data from data row and calculates statistical weight based on tmid error
        :param row:
        """
        self.epoch = row[0]
        self.tmid = row[1]
        self.tmid_err = row[2]
        try:
            self.weight = 1/self.tmid_err
        except ZeroDivisionError:
            self.weight = 0

    def __str__(self):
        return str(self.epoch)+' '+str(self.tmid)


def prop_to_date(row, end):
    """
    Calculate the propagated uncertainty based on the current data available
    :param row: array of data containing the information needed to propagate a timing uncertainty
    :param end: date to be propagated to as datetime object
    :return:
    """
    import numpy as np
    import julian

    # set date of ARIEL launch and convert to JD
    end_jd = julian.to_jd(end, fmt='jd') - 2400000
    # extract data
    period = float(row[1])
    period_err = float(row[2])
    tmid = row[3]
    tmid_err = row[4]
    duration = float(row[5])
    err_tot = 'NULL'
    current = tmid
    count = 0  # counter for epochs required

    if tmid_err is not None and period_err is not None:  # check that all information required is present
        if current < end_jd:  # check date is before launch
            while current < end_jd:  # loop towards it
                count += 1  # count epochs
                current += period  # increment date by period
            # propagate error by the number of epochs found and convert to percentage
            err_tot = np.sqrt(tmid_err * tmid_err + count * count * period_err * period_err)
            percent = err_tot * 24 * 60 / duration * 100
            # check for total loss
            if percent >= 100.0:
                loss = True
            else:
                loss = False
        else: # if date is after launch, just store current values
            print('End of sim reached')
            err_tot = tmid_err
            percent = err_tot*24*60/duration*100
            # check for total loss
            if percent >= 100.0:
                loss = True
            else:
                loss = False
    else: # if data missing, insert high values
        percent = 1000
        loss = True

    return err_tot, percent, loss


def prop_to_threshold(row, threshold, sim_end):
    import numpy as np
    import julian
    sim_end_jd = julian.to_jd(sim_end, fmt='jd') - 2400000
    # extract data
    period = float(row[1])
    period_err = float(row[2])
    tmid = row[3]
    tmid_err = row[4]
    err_tot = tmid_err
    count = 0  # counter for epochs required
    expiry_jd = 0
    days_to_threshold = 0
    if tmid_err is not None and period_err is not None:
        if err_tot < threshold/24/60:
            while err_tot < threshold/24/60:
                count += 1
                err_tot = np.sqrt(tmid_err * tmid_err + count * count * period_err * period_err)
                days_to_threshold = count*period
                if tmid + days_to_threshold > sim_end_jd:
                    expiry_jd = 100000
                    return expiry_jd
            expiry_jd = tmid + days_to_threshold
    print(row[0], row[-1], tmid, expiry_jd)
    return expiry_jd


