"""
This file contains the methods used to calculate various parameters needed by the PyBlastGraph program to produce
it's 2-D or 3-D graphs.
"""
import numpy as np


def restriction(text: str, a: float) -> float:
    """
    The restriction method calculates the restriction index of a blast site.

    :param text: str, the value of the corresponding radio button (indicating from which variable the restriction index
                                                                   will be calculated).
    :param a: float, the angle value in degrees.
    :return: float, the value of the restriction index.
    """
    f = None
    if text == "Angle":
        f = 3 / (3 + np.tan((np.pi / 2) - (a * np.pi / 180)))
    elif text == "Restriction Index":
        pass
    elif text == "Slope Ratio":
        f = 3 / (3 + a)

    return f


def angle_from_restriction(text: str, n: float) -> float:
    """
    The angle_from_restriction method calculates the slope of the face of a blast site.

    :param text: str, the value of the corresponding radio button (indicating from which variable the angle will be
                                                                   calculated).
    :param n: float, the face's ratio of inclination (i.e. 1/3 --> 0.3333).
    :return: float, the value of the face's angle in degrees.
    """
    a = None
    if text == "Angle":
        pass
    elif text == "Restriction Index":
        a = np.rad2deg((np.pi / 2) - np.arctan((3 / n) - 3))
    elif text == "Slope Ratio":
        a = np.rad2deg((np.pi / 2) - np.arctan(n))
    return a


def burden_langefors(d: float, p: float, s: float, c: float, f: float, sb: float) -> float:
    """
    The burden_langefors method calculates the burden that is to be applied in the blast's site calculations.

    :param d: float, the diameter of the drill hole (in mm or inches).
    :param p: float, the density of the explosive that is used.
    :param s: float, the relative strength of the explosive that is used.
    :param c: float, the rock constant of the blast site.
    :param f: float, the restriction index of the blast site.
    :param sb: float, the burden/spacing ratio.
    :return: float, the burden.
    """
    burden = (d / 33) * ((p * s) / (c * f * sb)) ** 0.5
    return burden


def face_length(k: float, a: float) -> float:
    """
    The face_length method calculates the length of the blast site face (in m or ft).

    :param k: float, the height of the face.
    :param a: float, the angle of the face's inclination.
    :return: float, the length of the face.
    """
    ab = k / np.cos((np.pi / 2) - (a * np.pi / 180))
    return ab


def error_collar(d: float) -> float:
    """
    The error_collar method calculates the deviation that will occur when drilling the holes and will be applied to the
    burden in order to correct it.

    :param d: float, the drill hole diameter
    :return: float, the error that will occur (in mm or inches).
    """
    e1 = d / 1000
    return e1


def error_drilling(h: float) -> float:
    """
    The error_drilling method calculates the deviation that will occur when drilling the holes and will be applied to
    the burden in order to correct it.

    :param h: float, the height of the drill hole
    :return: float, the error that will occur (in mm or inches).
    """
    e2 = 0.03 * h
    return e2


def practical_burden(burden: float, e1: float, e2: float) -> float:
    """
    The practical_burden method calculates the burden to be applied in practice after the corrects are applied to the
    theoretical value (from the burden_langefors method).

    :param burden: float, the theoretical burden from the burden_langefors method
    :param e1: float, the error that will occur (in mm or inches).
    :param e2: float, the error that will occur (in mm or inches).
    :return: float, the practical burden (in m or ft).
    """
    bp = burden - e1 - e2
    return bp


def bottom_length(kb: float, burden: float) -> float:
    """
    The bottom_length method calculates the length of the bottom explosive charge in a blast hole.

    :param kb: float, an index indicating the ratio of bottom to the whole length.
    :param burden: float, the theoretical burden.
    :return: float, the length of the bottom charge in a blast hole (in m or ft).
    """
    hb = (kb + 1) * burden
    return hb


def column_length(h: float, ho: float, hb: float) -> float:
    """
    The column_length method calculates the length of the column (or middle part) explosive charge in a blast  hole.

    :param h: float, the length of the whole drill hole.
    :param ho: float, the length of the stemming (top part with no explosives) of the drill hole.
    :param hb: float, the length of the bottom charge in a blast hole.
    :return: float, the length of the column charge.
    """
    hc = h - ho - hb
    return hc


def single_hole_volume(bp: float, spacing: float, k: float) -> float:
    """
    The single_hole_volume method calculates the volume of rock that a single blast hole will shutter.

    :param bp: float, the practical burden.
    :param spacing: float, the horizontal distance between adjacent drill holes.
    :param k: float, the height of the blast's site's face.
    :return: float the volume of rock that a single blast hole will shutter
    """
    vd = bp * spacing * k
    return vd


def rows_needed(l_max: float, bp: float) -> int:
    """
    The rows_needed method calculates the minimum amount of rows of drill holes needed in order to shutter the demanded
    rock volume.

    :param l_max: float, the maximum length of the blast site (from face to back wall).
    :param bp: float, the practical burden.
    :return: int, the minimum amount of rows of drill holes needed in this case.
    """
    rows = np.floor_divide(l_max, bp)
    return rows


def linear_density(e1: float, p: float) -> float:  # (Not used in PyBlastGraph)
    """
    The linear_density method calculates the linear density of the explosives load inside the drill hole.
    :param e1: float, the diameter of the drill hole.
    :param p: float the density of the explosive.
    :return: float, the linear density of the explosives.
    """
    lb = (np.pi * (e1 ** 2) * p) / 4
    return lb


def simple_sum(a, b):
    """
    A method to calculate a simple sum of two numbers.

    :param a:
    :param b:
    :return:
    """
    c = a + b
    return c


def simple_sub(a, b):
    """
    A method to calculate a simple subtraction of two numbers.

    :param a:
    :param b:
    :return:
    """
    c = a - b
    return c


def simple_multiplication(a, b):
    """
    A method to calculate a simple multiplication of two numbers.

    :param a:
    :param b:
    :return:
    """
    c = a * b
    return c


def simple_div(a, b):
    """
    A method to calculate a simple division of two numbers.

    :param a:
    :param b:
    :return:
    """
    c = a / b
    return c


def ceil_div(a: float, b: float) -> int:
    """
    The ceil_div method calculates the division of two numbers rounded to the next integer.

    :param a: float, the value of the numerator.
    :param b: float, the value of the denominator.
    :return: int, the result of the division rounded to the next integer.
    """
    c = np.ceil(a / b)
    return c


# # # # # # FRAGMENTATION AND RESTRICTIONS METHODS # # # # # # # # (Not used in PyBlastGraph)


def explosive_mass_index(qe: float, qa: float, ve: float, va: float) -> float:
    s = (5 / 6) * (qe / qa) + (1 / 6) * (ve / va)
    return s


def relative_explosive_power(s: float) -> float:
    e = s * 100
    return e


def relative_explosives_power(wqb: float, e: float, wc: float, ec: float) -> float:
    e_mean = (wqb * e + wc * ec) / (wqb + wc)
    return e_mean


def fragmentation(rf: float, vd: float, wqd: float, e: float) -> float:
    xm = rf * ((vd / wqd) ** 0.8) * (wqd ** (1 / 6)) * ((e / 115) ** (-19 / 30))
    return xm


def kuzram_optimization(rf: float, wqd: float, xm: float, e: float) -> float:
    q = (((rf * wqd ** (1 / 6)) / xm) * (e / 115) ** (-19 / 30)) ** (5 / 4)
    return q


def burden_kuzram(wqd: float, sb: float, q: float, k: float) -> float:
    burden = np.sqrt(wqd / (sb * q * k))
    return burden


def particle_velocity_restriction(kv: float, rr: float, w: float, nv: float) -> float:
    v = kv * (rr / (w ** (1 / 2))) ** (-nv)
    return v


def allowed_mass_v(rr: float, kv: float, v: float, nv: float) -> float:
    wapv = (rr * ((kv / v) ** (-1 / nv))) ** 2
    return wapv


def shockwave_restriction(kdb: float, rr: float, w: float, ndb: float) -> float:
    sp = kdb * (rr / (w ** (1 / 3))) ** (-ndb)
    return sp


def allowed_mass_sp(rr: float, kdb: float, sp: float, ndb: float) -> float:
    wadb = (rr * ((kdb / sp) ** (-1 / ndb))) ** 3
    return wadb


def decibel(sp: float, po: float) -> float:
    # po = 2.0 * (10 ** -9)
    db = 20 * np.log10(sp / po)
    return db


# # # # # # DELAY PATTERN METHODS # # # # # # # # (Not used in PyBlastGraph)


def delay_successively(nr):
    j = np.arange(1, nr + 1)
    a = j
    return a


def delay_per_row(r, c, stag, y):
    a = np.ones((r, c))  # The initial "Matrix" that will store the blast queue

    for i in range(1, r):
        a[i, :] = a[i-1, :] + 1

    a[:, 0] = a[:, 0] + 1
    a[:, c-1] = a[:, c-2] + 1

    if stag:
        for jj in range(r):
            if jj % 2 == 0:
                a[jj, c - 2] = a[jj, c - 3] + 1

        for jj in range(r):
            for i in range(c):
                if y[jj, i] == 0:
                    a[jj, i] = 0

    return a


def delay_echelon(c, r, stag, y):
    a = np.ones((r, c))  # The initial "Matrix" that will store the blast queue

    for i in range(1, c):
        a[0, i] = a[0, 0] + i

    for i in range(1, r):
        a[i, :] = a[0, :] + i

    if stag:
        for jj in range(r):
            for i in range(c):
                if y[jj, i] == 0:
                    a[jj, i] = 0

    return a


def delay_v_shape(r, c, stag, y):
    a = np.ones((r, c))          # The initial "Matrix" that will store the blast queue
    m = int(np.floor((c+1)/2))   # The closest possible column to the middle of the coordinates "matrices"
    a[:, m - 1] = np.arange(1, r + 1)

    for i in range(m - 2, -1, -1):
        a[:, i] = a[:, i + 1] + 1

    for i in range(m, c):
        a[:, i] = a[:, i - 1] + 1
    if stag:
        for jj in range(r):
            for i in range(c):
                if y[jj, i] == 0:
                    a[jj, i] = 0

        a[:, -1] -= 1

        for jj in range(r):
            for i in range(c):
                if a[jj, i] < 0:
                    a[jj, i] = 0

    return a


def delay_double_v(r, c, stag, y):
    a = np.ones((r, c))          # The initial "Matrix" that will store the blast queue
    m = int(np.floor((c+1)/2))   # The closest possible column to the middle of the coordinates "matrices"
    m1 = int(np.floor(m/2))      # The middle column between column '0' and column 'm' as close as possible
    m2 = int(np.floor((m+c)/2))  # The middle column between column 'm' and the last column as close as possible

    a[:, m1-1] = np.arange(1, r+1)
    a[:, m2-1] = np.arange(1, r+1)

    for i in range(m1-2, -1, -1):
        a[:, i] = a[:, i+1] + 1

    for i in range(m1, m):
        a[:, i] = a[:, i-1] + 1

    for i in range(m2-2, m-1, -1):
        a[:, i] = a[:, i+1] + 1

    for i in range(m2, c):
        a[:, i] = a[:, i-1] + 1

    if stag:
        for jj in range(r):
            for i in range(c):
                if y[jj, i] == 0:
                    a[jj, i] = 0

        for jj in range(r):
            if a[jj, -1] != 0:
                a[jj, -1] = a[jj, -2]

    return a


def delay_double_int_v(r, c, stag, y):
    a = np.ones((r, c))          # The initial "Matrix" that will store the blast queue
    m = int(np.floor((c+1)/2))   # The closest possible column to the middle of the coordinates "matrices"
    m1 = int(np.floor(m/2))      # The middle column between column '0' and column 'm' as close as possible
    m2 = int(np.floor((m+c)/2))  # The middle column between column 'm' and the last column as close as possible

    a[:, m1 - 1] = np.arange(1, r + 1)

    for i in range(m1 - 2, -1, -1):
        a[:, i] = a[:, i + 1] + 1

    for i in range(m1, m):
        a[:, i] = a[:, i - 1] + 1

    last = np.max(a) + 1
    a[:, m2 - 1] = np.arange(last, r + last)

    for i in range(m2 - 2, m - 1, -1):
        a[:, i] = a[:, i + 1] + 1

    for i in range(m2, c):
        a[:, i] = a[:, i - 1] + 1

    if stag:
        for jj in range(r):
            for i in range(c):
                if y[jj, i] == 0:
                    a[jj, i] = 0

        for jj in range(r):
            if a[jj, -1] != 0:
                a[jj, -1] = a[jj, -2]

    return a


def queue_table(x, y, a, r, c, stag, pattern):
    if not stag:
        # reshape "matrices" into arrays
        x = np.reshape(x, r*c)
        y = np.reshape(y, r*c)
        queue = np.reshape(a, r*c)  # Transform tha 'a' matrix into the 'queue' array (for understanding purposes)

        total1 = np.lexsort((y, x, queue))  # Sort the X, Y coordinates by the queue

        # Create the final "matrix" containing all data
        total1 = np.array([[x[i], y[i], queue[i].astype(int)]for i in total1])
    else:
        # reshape "matrices" into arrays
        row, col = np.shape(x)
        x = np.reshape(x, row*col)
        y = np.reshape(y, row*col)

        # If delay pattern is the first one only
        queue = a.astype(int)

        # If delay pattern is anything else
        if pattern != "Successively":
            queue = np.reshape(a, r*c)  # Transform tha 'a' matrix into the 'queue' array (not actually needed)
            queue = np.delete(queue, np.where(y == 0))  # Delete the zero entries

        x = np.delete(x, np.where(y == 0))  # Delete the zero entries
        y = np.delete(y, np.where(y == 0))  # Delete the zero entries
        print(f"x = {np.shape(x)}, y = {np.shape(y)}, queue = {np.shape(queue)}")

        total1 = np.lexsort((x, y, queue))  # Sort the X, Y coordinates by the queue

        # Create the final "matrix" containing all data
        total1 = np.array([[x[i], y[i], queue[i].astype(int)]for i in total1])

    return x, y, total1
