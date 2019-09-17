import math
from PyAstronomy import pyasl
import datetime
import numpy as np


class Error(Exception):
    pass


class NeverVisibleError(Error):
    pass


class AlwaysVisibleError(Error):
    pass

def get_gmst(j_date):

    # input: julian date
    # return GMST for a given julian date

    du = j_date - 2451545
    du_mod = math.fmod(du, 1)
    T = (j_date - 2451545) / 36525
    gmst_sec = 86400 * (
                0.7790572732640 + 0.00273781191135448 * du + du_mod) + 0.00096707 + 307.47710227 * T + 0.092772113 * (
                           T ** 2) - 0.0000000293 * (T ** 3) + 0.00000199707 * (T ** 4) - 0.000000002453 * (T ** 5)
    gmst_deg = math.fmod(gmst_sec / 3600, 24)
    return gmst_deg


def get_gmst_utc_diff(date):

    # input: UTC datetime
    # get gmst/utc difference for a given date
    # you can then get lst by subtracting EAST longitude

    jd = pyasl.jdcnv(date)  # find Julian date for given date in UTC
    gmst_dec = get_gmst(jd)  # get gmst in decimal format
    date_dec = date.hour + date.minute / 60 + date.second / 3600 # get UTC for date in decimal format
    gmst_utc_diff = gmst_dec - date_dec # get gmst - utc difference
    if gmst_utc_diff < 0:
        gmst_utc_diff += 24
    return gmst_utc_diff


def get_set_rise(start_date, ra, dec, lon, lat, h=0, skip_negative_check=False):
    import math
    import datetime

    factor = (math.sin(h * math.pi / 180) - math.sin(dec * math.pi / 180) * math.sin(
        lat * math.pi / 180)) / (math.cos(dec * math.pi / 180) * math.cos(lat * math.pi / 180))

    if not factor < -1 and not factor > 1:

        HA = math.acos(factor) * 180 / math.pi

        rise = ra - HA / 15. + 24 # shift rise by 24h, so we have set first, then rise
        set = ra + HA / 15.

        if set < 0:
            rise += 24
            set += 24

        # adjust using gmst / utc difference for start date (midnight, need further adjustment - see later)
        gmst_utc_diff = get_gmst_utc_diff(start_date)  # this is only approximate as it's at midnight of start date
        rise_lst = rise - gmst_utc_diff
        set_lst = set - gmst_utc_diff
        if set_lst < 0:
            rise_lst += 24
            set_lst += 24

        # refine lst adjustment for rise and set times (get gmst/utc difference at exactly set, rise)
        rise_lst_adj_tmp = 0
        set_lst_adj_tmp = 0
        gmst_utc_diff_rise = gmst_utc_diff
        gmst_utc_diff_set = gmst_utc_diff
        diff = 99
        while diff > 1: # until accuracy is less than a second
            set_final_tmp = set - gmst_utc_diff_set - lon / 15
            rise_final_tmp = rise - gmst_utc_diff_rise - lon / 15
            if set_final_tmp < 0:
                set_final_tmp += 24
                rise_final_tmp += 24
            date_set_UTC = start_date + datetime.timedelta(seconds=set_final_tmp * 60 * 60)
            date_rise_UTC = start_date + datetime.timedelta(seconds=rise_final_tmp * 60 * 60)
            gmst_utc_diff_rise = get_gmst_utc_diff(date_rise_UTC)
            gmst_utc_diff_set = get_gmst_utc_diff(date_set_UTC)
            rise_lst_adj = rise - gmst_utc_diff_rise
            set_lst_adj = set - gmst_utc_diff_set
            if set_lst_adj < 0:
                rise_lst_adj += 24
                set_lst_adj += 24
            diff = np.abs(rise_lst_adj_tmp - rise_lst_adj) + np.abs(set_lst_adj_tmp - set_lst_adj)
            rise_lst_adj_tmp = rise_lst_adj
            set_lst_adj_tmp = set_lst_adj
        rise_lst = rise_lst_adj
        set_lst = set_lst_adj

        # adjust for longitude
        rise_final = rise_lst - lon / 15
        set_final = set_lst - lon / 15

        if set_final < 0 and not skip_negative_check:

            # set is the previous day. Thefore get set and rise times for next day, then add 24h
            start_date_plus1day = start_date + datetime.timedelta(days=1)
            set_plus1day, rise_plus1day = get_set_rise(start_date_plus1day, ra, dec, lon, lat, h, skip_negative_check=True)
            return set_plus1day + 24, rise_plus1day + 24

        else:
            return set_final, rise_final

    else:

        if factor < 1:
            # target is always observable
            raise Warning('AlwaysVisible')

        elif factor > 1:
            # target is never observable
            raise Warning('NeverVisible')


def get_rise_set(start_date, ra, dec, lon, lat, h=0, skip_negative_check=False):
    factor = (math.sin(h * math.pi / 180) - math.sin(dec * math.pi / 180) * math.sin(
        lat * math.pi / 180)) / (math.cos(dec * math.pi / 180) * math.cos(lat * math.pi / 180))

    if not factor < -1 and not factor > 1:

        HA = math.acos(factor) * 180 / math.pi

        rise = ra - HA / 15.
        set = ra + HA / 15.

        if rise < 0:
            rise += 24
            set += 24

        # adjust using gmst / utc difference for start date (midnight, need further adjustment - see later)
        gmst_utc_diff = get_gmst_utc_diff(start_date) # this is only approximate
        rise_lst = rise - gmst_utc_diff
        set_lst = set - gmst_utc_diff
        if rise_lst < 0:
            rise_lst += 24
            set_lst += 24

        # refine lst adjustment for rise and set times (get gmst/utc difference at exactly set, rise)
        rise_lst_adj_tmp = 0
        set_lst_adj_tmp = 0
        gmst_utc_diff_rise = gmst_utc_diff
        gmst_utc_diff_set = gmst_utc_diff
        diff = 99
        while diff > 1: # until accuracy is less than a second
            set_final_tmp = set - gmst_utc_diff_set - lon / 15
            rise_final_tmp = rise - gmst_utc_diff_rise - lon / 15
            if rise_final_tmp < 0:
                set_final_tmp += 24
                rise_final_tmp += 24
            date_rise_UTC = start_date + datetime.timedelta(seconds = rise_final_tmp * 60 * 60)
            date_set_UTC = start_date + datetime.timedelta(seconds = set_final_tmp * 60 * 60)
            gmst_utc_diff_rise = get_gmst_utc_diff(date_rise_UTC) # get gmst-utc
            gmst_utc_diff_set = get_gmst_utc_diff(date_set_UTC) # get gmst-utc
            rise_lst_adj = rise - gmst_utc_diff_rise
            set_lst_adj = set - gmst_utc_diff_set
            if rise_lst_adj < 0:
                rise_lst_adj += 24
                set_lst_adj += 24
            diff = np.abs(rise_lst_adj_tmp-rise_lst_adj) + np.abs(set_lst_adj_tmp-set_lst_adj)
            rise_lst_adj_tmp = rise_lst_adj
            set_lst_adj_tmp = set_lst_adj
        rise_lst = rise_lst_adj
        set_lst = set_lst_adj

        # adjust for longitude
        rise_final = rise_lst - lon / 15
        set_final = set_lst - lon / 15

        if rise_final < 0 and not skip_negative_check:

            # rise is the previous day. Thefore get rise and set times for next day, then add 24h
            start_date_plus1day = start_date + datetime.timedelta(days = 1)
            rise_plus1day, set_plus1day = get_rise_set(start_date_plus1day, ra, dec, lon, lat, h, skip_negative_check=True)
            return rise_plus1day+24, set_plus1day+24

        else:
            return rise_final, set_final

    else:

        if factor < 1:
            # target is always observable
            raise AlwaysVisibleError

        elif factor > 1:
            # target is never observable
            raise NeverVisibleError


def sun_set_rise(start_date, lon, lat, sundown):
    from PyAstronomy import pyasl

    jd_set = pyasl.jdcnv(start_date)  # find Julian date
    sunpos = pyasl.sunpos(jd_set)
    sun_ra, sun_dec = sunpos[1][0], sunpos[2][0]

    try:
        set, rise = get_set_rise(start_date, sun_ra/15, sun_dec, lon, lat, h=sundown)
        set_date_final = start_date + datetime.timedelta(seconds=set * 60 * 60)
        rise_date_final = start_date + datetime.timedelta(seconds=rise * 60 * 60)

        # TODO we should actually adjust for RA and DEC at rise / set

        return set_date_final, rise_date_final

    except:

        pass


def target_rise_set(date, ra, dec, lon, lat, mintargetalt):

    rise, set = get_rise_set(date, ra / 15, dec, lon, lat, h=mintargetalt)

    rise_dt = date + datetime.timedelta(seconds=rise * 60 * 60)
    set_dt = date + datetime.timedelta(seconds=set * 60 * 60)

    return rise_dt, set_dt
