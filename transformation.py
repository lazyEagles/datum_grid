
import re
import math

BLOCK_X = {
    90000: "J M P".split(),
    170000: "A D G K N Q T V".split(),
    250000: "B E H L O R U W".split(),
    330000: "C F I".split()
}

BLOCK_Y = {
    2750000: "A B C".split(),
    2700000: "D E F".split(),
    2650000: "G H I".split(),
    2600000: "J K L".split(),
    2550000: "M N O".split(),
    2500000: "P Q R".split(),
    2450000: "T U".split(),
    2400000: "V W".split()
}

DECODE_G_X = {
    y: x for x, ys in BLOCK_X.items() for y in ys
}

DECODE_G_Y = {
    y: x for x, ys in BLOCK_Y.items() for y in ys
}

DECODE_R = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5, "G": 6, "H": 7}
DECODE_S = {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4}


def to_twd67tm2_from(taipower_grid):
    """
    transform taipower_grid to twd67tm2

    ref:
    http://www.sunriver.com.tw/grid_taipower.htm
    https://wiki.osgeo.org/wiki/Taiwan_Power_Company_grid
    """

    p = re.compile(r"([a-zA-Z])(\d\d)(\d\d) *([a-zA-Z])([a-zA-z])(\d)(\d)((\d)(\d))?")
    m = p.search(taipower_grid)

    G_x = DECODE_G_X[m.group(1)]
    G_y = DECODE_G_Y[m.group(1)]
    PP = int(m.group(2))
    QQ = int(m.group(3))
    R = DECODE_R[m.group(4)]
    S = DECODE_S[m.group(5)]
    T = int(m.group(6))
    U = int(m.group(7))
    V = int(m.group(9)) if m.group(9) else 0
    M = int(m.group(10)) if m.group(10) else 0

    x = G_x + PP * 800 + R * 100 + T * 10 + V
    y = G_y + QQ * 500 + S * 100 + U * 10 + M

    return x, y

def to_twd67latlon_from(twd67tm2):
    """
    tranform twd67tm2 to twd67latlon

    ref:
    https://en.wikipedia.org/wiki/Transverse_Mercator_projection
    http://www.sunriver.com.tw/grid_tm2.htm
    http://www.ihsenergy.com/epsg/guid7.pdf
    """
    E, N = twd67tm2

    a = 6378160
    f = 1 / 298.25
    e_square = 2*f - f**2
    e1 = (1-math.sqrt(1-e_square)) / (1+math.sqrt(1-e_square))
    k0 = 0.9999

    FE = 250000
    FN = 0

    M0 = 0
    M1 = M0 + (N - FN)/k0
    mu1 = M1 / (a*(1- e_square/4 - 3*e_square**2/64 - 5*e_square**3/256))
    phi1 = (
        mu1
        + (3*e1/2 - 27*e1**3/32)*math.sin(2*mu1)
        + (21*e1**2/16 - 55*e1**4/32)*math.sin(4*mu1)
        + (151*e1**3/96)*math.sin(6*mu1)
        + (1097*e1**4/512)*math.sin(8*mu1)
    )
    v1 = a / math.sqrt(1-e_square*math.sin(phi1)**2)
    rho1 = a * (1-e_square) / (1-e_square*math.sin(phi1)**2)**1.5
    T1 = math.tan(phi1)**2
    e_p_square = e_square / (1-e_square)
    C1 = e_p_square * math.cos(phi1)**2
    D = (E-FE) / (v1*k0)

    phi = (
        phi1
        - (v1*math.tan(phi1)/rho1)*(
            D**2/2
            - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*e_p_square)*D**4/24
            + (
                61 + 90*T1 + 298*C1 + 45*T1**2 - 252*e_p_square - 3*C1**2
            )*D**6/720
        )
    )

    lamb0 = 121 * math.pi / 180

    lamb = (
        lamb0
        + (
            D
            - (1+2*T1+C1)*D**3/6
            + (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*e_p_square + 24*T1**2)*D**5/120
        )/math.cos(phi1)
    )

    lat = phi * 180 / math.pi
    lon = lamb * 180 / math.pi

    return lat, lon


def main():
    tg1 = "G8150HD7812"
    tg2 = "B8146CC58"
    print(to_twd67latlon_from(to_twd67tm2_from(tg1)))
    print(to_twd67latlon_from(to_twd67tm2_from(tg2)))

if __name__ == "__main__":
    main()