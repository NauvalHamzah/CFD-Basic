import math


# --------------------------------------------------------------------------------- #

def NACA00TT(x, t, c):
    yt = 5*t*c*(0.2969*math.pow(x,0.5) - 0.126*x - 0.3516*math.pow(x,2) + 0.2843*math.pow(x,3) - 0.1036*math.pow(x,4))

    return yt

def NACAMPTT(x, m, p, t, c):
    yt = 5*t*(0.2969*math.pow(x,0.5) - 0.126*x - 0.3516*math.pow(x,2) + 0.2843*math.pow(x,3) - 0.1036*math.pow(x,4))
    
    if x < p:
        yc = m*(2*p*x - math.pow(x,2))/math.pow(p,2)
    else:
        yc = m*(1 - 2*p + 2*p*x - math.pow(x,2))/math.pow((1-p),2)
    
    if x < p:
        gyc = 2*m*(p - x)/math.pow(p,2)
    else:
        gyc = 2*m*(p - x)/math.pow((1-p),2)
    teta = math.atan(gyc)

    xu = c*(x - yt*math.sin(teta))
    xl = c*(x + yt*math.sin(teta))
    yu = c*(yc + yt*math.cos(teta))
    yl = c*(yc - yt*math.cos(teta))

    return xu, xl, yu, yl

# --------------------------------------------------------------------------------- #