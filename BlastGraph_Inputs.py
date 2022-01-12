"""
This is an input example for the development of the 2D and 3D graphs for PyBlastGraph.

Returns the inputs needed in PyBlastGraph (burden, spacing, rows ect).

"""
import BlastGraph_Equations as eQ
import numpy as np


print("Entered DesignParameters")

# Inputs
a = 72      # angle value in degrees
d = 76      # hole diameter in mm
p = 0.8     # explosive's density in g/cm^3
s = 0.84    # explosive's relative strength (dimensionless)
c = 0.45    # the rock constant (dimensionless)
sb = 1.25   # the spacing to burden ratio (dimensionless)
k = 10      # the height of the site's face in m
lmax = 30   # the maximum length of the blast site in m
vp = 10000  # the demanded rock volume to be shuttered in m^3
kb = 0.3    # the bottom charge ratio (dimensionless)
kt = 1      # the stemming ratio (dimensionless)
a1 = 70     # the angle between the two faces of the site (if there are two faces)

# Outputs
f = eQ.restriction("Angle", a)
b = eQ.burden_langefors(d, p, s, c, f, sb)
ab = eQ.face_length(k, a)
u = kb * b                                      # the subdrill length in m
h = ab + u                                      # the drill hole length in m
e1 = eQ.error_collar(d)
e2 = eQ.error_drilling(h)
bp = eQ.practical_burden(b, e1, e2)
sp = eQ.simple_multiplication(sb, bp)
ho = eQ.simple_multiplication(kt, b)
hb = eQ.bottom_length(kb, b)
hc = eQ.column_length(h, ho, hb)
vd = eQ.single_hole_volume(bp, sp, k)
nsimple = eQ.ceil_div(vp, vd)                   # the minimum amount of blast holes needed to produce the given volume
rows = eQ.rows_needed(lmax, bp)
perrow = eQ.ceil_div(nsimple, rows)             # the blast holes in each row
ntot = eQ.simple_multiplication(rows, perrow)   # the actual amount of blast holes
hcheck = ho + hc + hb                           # a check of the lengths (should be equal to h)

# PRESPLIT
presplitting = "y"
# presplitting = "n"
if presplitting == "y":
    # INPUTS
    dholeps = 35    # the diameter of the pre splitting holes in mm
    dexps = 28      # the cartridge diameter in mm
    pps = 1.2       # the density of the explosive in g/cm^3
    tps = 15.1      # the tensile strength of the rock in MPa
    vps = 4500      # the velocity of the explosion in m/s
    ktps = 20       # the stemming index

    # OUTPUT
    hpstot = k / np.cos((np.pi / 2) - (a * np.pi / 180))            # the length of the pre splitting drill hole in m
    dsexps = dexps / 1000                                           # transform the diameter from mm to m
    hops = ktps * dsexps                                            # the stemming length in m
    hps = hpstot - hops                                             # the length of the explosive column

    prps = ((0.000228 * pps * (vps ** 2)) / (1 + 0.8 * pps)) * \
           ((np.sqrt(k / hps)) * (dholeps / dexps)) ** 2.4          # the explosion pressure in Mpa
    dexps = dexps / 1000                                            # transform the diameter from mm to m
    sps = dexps * ((prps / (tps * 5)) + 1)                          # the spacing of the drill holes in m
    rl = rows * bp                                                  # the actual length of the site in m
    widthtrue = perrow * sp                                         # the actual width of the site in m
    nps = np.ceil(rl / sps)                                         # the drill holes needed for one side of the site
    npsback = nps                                                   # the drill holes needed for the other side

    nrps = npsback + nps                                            # total drill holes needed

    # The below are for the 2 faces case
    ltrue2 = bp * rows + (sp / 2)
    widthtrue2 = bp + (perrow * sp) - (sp / 2)
    rl2 = ltrue2 / np.cos(np.pi / 2 - a1 * np.pi / 180)
    nps2 = np.ceil(rl2 / sps)
    npsback2 = np.floor(widthtrue2 / sps)
    nrps2 = npsback2 + nps2

print("Exited DesignParameters")
