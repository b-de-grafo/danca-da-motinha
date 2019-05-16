import numpy as np

def runge_kutta(h,y_seta,y_zero, intervalo):
    y = [y_zero]
    x = list(np.arange(intervalo[0], intervalo[1], h))
    for k in range(len(x)-1):
        y.append(None)
        y[k+1] = y[k] + (h/2) * (y_seta(x[k], y[k]) + y_seta(x[k] + h, y[k]+h*y_seta(x[k], y[k])))
    return(y)
    #y_um = y_zero + (h/2) * (y_seta(x_zero,y_zero) + y_seta(x_zero+h, y_zero+h*y_seta(x_zero,y_zero))

def y_seta_ex(x,y):
    return -y+x-2

print(runge_kutta(0.1,y_seta_ex,2,[0,0.5]))

