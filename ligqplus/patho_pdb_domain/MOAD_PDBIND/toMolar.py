import math

def toMolar(valor,unidad):
    """Toma una concetracion y su unidad para convertirla a Molar y calcular su logaritmo negativo"""
    units_convertion={
        "M^-1":0.1,
        "M":1,
        "mM":0.001,
        "uM":0.000001,
        "nM":0.000000001,
        "pM":0.000000000001,
        "fM":0.000000000000001
    }

    return round(math.log(valor*units_convertion[unidad])*-1,3)
