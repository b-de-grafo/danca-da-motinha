import numpy as np
import matplotlib.pyplot as plt


def runge_kutta(h, y_seta, y_zero, intervalo):
    y = [y_zero]
    x = list(np.arange(intervalo[0], intervalo[1], h))
    for k in range(len(x)-1):
        y.append(None)
        y[k+1] = y[k] + (h/2) * (y_seta(x[k], y[k]) + y_seta(x[k] + h, y[k]+h*y_seta(x[k], y[k])))
    return x, y


def passo_multiplo(h, y_seta, y_zero, intervalo, inicializacao):
    y = [y_zero]
    x = list(np.arange(intervalo[0], intervalo[1], h))
    
    y.append(inicializacao(h, y_seta, y_zero, [intervalo[0], intervalo[0]+2*h])[1][1])    
    
    for k in range(1, len(x)-1):
        y.append(None)
        y[k+1] = y[k] + (h/2) * (3 * y_seta(x[k], y[k]) - y_seta(x[k-1], y[k-1]))
    
    return x, y

def y_seta_ex(x, y):
    return -y+x-2

x1, y1 = runge_kutta(0.1, y_seta_ex, 2, [0, 5])
fig, axs = plt.subplots(1, 1)
axs.plot(x1, y1, color = "red")

x2, y2 = passo_multiplo(0.1, y_seta_ex, 2, [0, 5], runge_kutta)
axs.plot(x2, y2, color = "blue")

plt.show()