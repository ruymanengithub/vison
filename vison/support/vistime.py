#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""

Accesory library: time related operations

Created on Tue Oct 10 15:08:28 2017

:author: Ruyman Azzollini

"""

# IMPORT STUFF
import datetime
# END IMPORT

dtobj_default = datetime.datetime(1980, 2, 21, 7, 0, 0)  # early riser


def get_time_tag():
    """ """
    t = datetime.datetime.now()
    s = t.strftime('%Y%m%d_%H%M%S')
    return s


def get_dtobj(DT):
    """ """

    date = DT[0:DT.index('D')]
    y2d = int(date[4:6])
    if y2d < 20:
        century = 2000
    else:
        century = 1900
    dd, MM, yy = int(date[0:2]), int(date[2:4]), y2d + century

    time = DT[DT.index('D') + 1:-1]

    hh, mm, ss = int(time[0:2]), int(time[2:4]), int(time[4:6])

    dtobj = datetime.datetime(yy, MM, dd, hh, mm, ss)

    return dtobj
