# coding: utf-8
def Kcorr(ty,z):
    if 0.4 < z < 1.0:
        if ty == 1:
            return 1.2*z**2
        elif ty == 2: 
            return 0.9*z**3
        elif ty == 3:
            return 0.5*z**3
        else:
            return -0.25
    elif 0.0 < z < 0.4:
        if ty == 1:
            return (1.2*0.4**2)/0.4 * z
        elif ty == 2:
            return (0.9*0.4**3)/0.4 * z
        elif ty == 3:
            return (0.5*0.4**3)/0.4 * z
        else:
            return -0.25/0.4 * z
    else: return np.nan
        
