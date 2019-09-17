#################################################################
# Code to query databases for observation and target data       #
# Called by methods within database_generator.py                #
#                                                               #
# Hamish Caines 07-2019                                         #
#################################################################


class Target:
    """
    Target object, contains empty variables for information obtained from queries
    Functions within make the queries and populate variables
    """
    def __init__(self, name):
        """
        Constructor for new information about a single target
        :param name: Name of target planet
        """
        self.name = name
        self.star = name[:-1]  # extract star name
        self.planet = name[-1]  # extract planet id
        self.ETD_obs = []
        self.query_depth = None
        self.query_depth_err = None
        self.query_lastob = None
        self.query_lastob_err = None
        self.query_last_epoch = None

        self.depth_source = ''

        self.EXO_used = False

    def __str__(self):
        return self.star+self.planet

    def ETD_query(self):
        """
        Query the Exoplanet Target Database for observation data and store in object
        """
        import re
        import requests
        import numpy as np
        # build url for target
        url_base = 'http://var2.astro.cz/ETD/etd.php?'
        url = url_base + 'STARNAME=' + self.star + '&PLANET=' + self.planet

        web_html = str(requests.get(url).content)  # obtain html from page for target

        # obtain data from html using regular expressions
        obs_data = re.findall(
            "<tr valign=\\\\\\'top\\\\\\'><td>(\d+)<\/td><td class=\\\\\\'right\\\\\\'><b>([\d\.]+)<br\/><\/b>"
            "([\d\s\.\+\/-]+)<\/td><td>(\d+)<\/td><td>([+|\-\d\.\s]+)<\/td><td>([\d\.\s]+)\+\/-([\d\.\s]+)<\/td><td>"
            "([\d\.]+) ([\d\s\.\+\/-]+)<\/td><td>([\w]+)<\/td><td><b><a href=\\\\\\'etd-data\.php\?id=[\d]+\\\\\\' "
            "target=\\\\\\'[\w]+\\\\\\' title=\\\\\\'get data\\\\\\'>(\d)", web_html)

        # initialise depth counters
        depth_count = 0
        depth_total = 0
        depth_err_total = 0
        # initialise latest observation variables
        max_tmid = 0
        max_tmid_err = 0
        max_epoch = 0
        # iterate through observations
        for obs in obs_data:
            new_ob = Observation()
            new_ob.gen_from_ETD(obs)

            # check quality assigned by ETD
            if new_ob.quality >= 3:
                self.ETD_obs.append(new_ob)

                # extract depth data and add to counters
                if new_ob.depth_err is not None:
                    depth_total += new_ob.depth
                    depth_err_total += new_ob.depth_err * new_ob.depth_err
                    depth_count += 1

                # extract and store tmid data
                if new_ob.tmid > max_tmid:
                    max_tmid = new_ob.tmid
                    max_tmid_err = new_ob.tmid_err
                    max_epoch = new_ob.epoch

        # store tmid data in object
        if max_tmid != 0:
            self.query_lastob = max_tmid
            self.query_lastob_err = max_tmid_err
            self.query_last_epoch = max_epoch

        # calculate and store depth data
        try:
            self.query_depth = depth_total / depth_count
            self.query_depth_err = np.sqrt(depth_err_total) / depth_count
            self.depth_source = 'ETD'
        except ZeroDivisionError:
            pass

    def EXO_query(self):
        """
        Query exoplanets.org for missing data and store in object
        """
        import requests
        import re

        # build url for specific target, many cases
        url_base = 'http://exoplanets.org/detail/'
        if self.star[0].isdigit():
            number = True
            count = 0
            while number:
                count += 1
                number = self.star[count].isdigit()
            url = url_base + self.star[:count] + '_' + self.star[count:] + '_' + self.planet
        elif self.star[:4] == 'EPIC':
            url = url_base + self.star[:4] + '_' + self.star[4:] + '_' + self.planet
        elif 'GJ' in self.star:
            url = url_base + 'GJ' + '_' + self.star[2:] + '_' + self.planet
        elif 'HD' in self.star:
            url = url_base + 'HD' + '_' + self.star[2:] + '_' + self.planet
        elif 'LHS' in self.star:
            url = url_base + 'LHS' + '_' + self.star[3:] + '_' + self.planet
        elif 'NGTS' in self.star:
            url = url_base + 'NGTS' + '-' + self.star[4:] + '_' + self.planet
        elif 'PH' in self.star:
            url = url_base + 'PH-' + self.star[2:] + '_' + self.planet
        elif 'KPS' in self.star:
            url = url_base + 'KPS' + '-' + self.star[3:] + '_' + self.planet
        elif 'HIP' in self.star:
            url = url_base + 'HIP' + '_' + self.star[3:] + '_' + self.planet
        else:
            url = url_base + self.star + '_' + self.planet

        # obtain html from url
        web_html = str(requests.get(url).content)

        # check for valid planet
        invalid_check = re.findall("(Invalid Planet)", web_html)
        if invalid_check:
            print(url + ' is invalid')
        if not invalid_check:
            # read data
            data = dict(re.findall("\"([\w]+)\":\[([\d\w\.\"\s\/\:\;\%]+)\]", web_html))

            # extract required information and store in object
            try:
                if self.query_depth is None:
                    if data['DEPTH'] != 'null':
                        self.query_depth = float(data['DEPTH']) * 1000  # convert from unitless to mmag
                        self.depth_source = 'EXO'
                        if data['DEPTHUPPER'] != 'null':
                            self.query_depth_err = float(data['DEPTHUPPER']) * 1000  # convert from unitless to mmag

                if self.query_last_epoch is None:
                    if data['TT'] != 'null':
                        self.query_lastob = float(data['TT']) - 2400000
                        self.query_lastob_err = data['TTUPPER']
                        self.query_last_epoch = 0
            except KeyError:
                pass

    def write_query_data(self, cursor, db):
        """
        Store data from queries in database
        :param cursor: Cursor connected to database: cursor
        :param db: Database: database
        """
        import sqlite3
        # check for new data and store in database
        if self.query_depth is not None:
            cursor.execute('UPDATE TARGET_DATA SET Depth = '+str(self.query_depth)+' WHERE Name = \''+self.name+'\'')
            cursor.execute('UPDATE TARGET_DATA SET DepthSource = \''+self.depth_source+'\' WHERE Name = \''+self.name+'\'')
            if self.query_depth_err is not None:
                cursor.execute('UPDATE TARGET_DATA SET DepthErr = '+str(self.query_depth_err) + \
                         ' WHERE Name = \''+self.name+'\'')
        if self.query_last_epoch is not None:
            cursor.execute('UPDATE TARGET_DATA SET LastObs = '+str(self.query_lastob)+' WHERE Name = \''+self.name+'\'')
            cursor.execute(
                'UPDATE TARGET_DATA SET LastObsErr = ' + str(self.query_lastob_err) + ' WHERE Name = \'' + self.name + '\'')
            cursor.execute(
                'UPDATE TARGET_DATA SET LastEpoch = ' + str(self.query_last_epoch) + ' WHERE Name = \'' + self.name + '\'')

        # check number of observations and store in database
        if len(self.ETD_obs) != 0:
            cursor.execute('UPDATE TARGET_DATA SET NoOfObs = '+str(len(self.ETD_obs))+' WHERE Name = \''+self.name+'\'')
            # store each observation in table for target
            for ob in self.ETD_obs:
                try:
                    cursor.execute('INSERT INTO \''+self.name+'\' VALUES ('+ob.obnumber+', '+str(ob.epoch)+', '+str(ob.tmid)+\
                    ', '+str(ob.tmid_err)+', NULL, NULL, '+str(ob.depth)+', '+str(ob.depth_err)+', '+str(ob.duration)+', '+\
                    str(ob.duration_err)+', \'ETD\')')
                except sqlite3.IntegrityError:
                    print('Obs No. '+str(ob.obnumber)+' already exists')
        db.commit()


class Observation:
    """
    Observation object, contains empty variables for observation data from queries
    Populated by function within
    """
    def __init__(self):
        """
        Null constructor
        """
        self.obnumber = 0
        self.epoch = 0
        self.tmid = 0
        self.tmid_err = 0
        self.depth = 0
        self.depth_err = 0
        self.duration = 0
        self.duration_err = 0
        self.quality = 0
        self.OminusC = 0

    def gen_from_ETD(self, match):
        """
        Populute Observation object with data obtained for a single observation from ETD
        :param match: list of variables for a single observations
        """
        self.obnumber = match[0]
        self.tmid = float(match[1])  # HJD
        try:
            self.tmid_err = float(match[2].split()[1])  # HJD
        except ValueError:
            self.tmid_err = 0
        except IndexError:
            self.tmid_err = 0
        self.epoch = int(match[3])
        self.OminusC = float(match[4])  # days
        try:
            self.duration = float(match[5])  # mins
        except ValueError:
            pass
        try:
            self.duration_err = float(match[6].split()[0])  # mins
        except IndexError:
            pass
        self.depth = float(match[7])  # mmag
        try:
            self.depth_err = float(match[8].split()[1])  # mmag
        except IndexError:
            pass
        try:
            self.quality = int(match[10])
        except IndexError:
            pass
