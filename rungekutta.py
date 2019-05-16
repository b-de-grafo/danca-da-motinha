import numpy as np
import matplotlib.pyplot as plt


def runge_kutta(h, y_seta, y_zero, intervalo):
    y = [y_zero]
    x = list(np.arange(intervalo[0], intervalo[1], h))
    for k in range(len(x)-1):
        y.append(None)
        y[k+1] = y[k] + (h/2) * (y_seta(x[k], y[k]) + y_seta(x[k] + h, y[k]+h*y_seta(x[k], y[k])))
    return x, y

    # y_um = y_zero + (h/2) * (y_seta(x_zero,y_zero) + y_seta(x_zero+h, y_zero+h*y_seta(x_zero,y_zero))


def y_seta_ex(x, y):
    return -y+x-2

def grafico(x, y):
    vetor = []
    for i in range(len(x)):
        vetor.append([x[i], y[i]])

    print(vetor)

    fig, axs = plt.subplots(1, 1)
    axs.plot(vetor)
    plt.show()



x, y = runge_kutta(0.5, y_seta_ex, 2, [0, 5])

fig, axs = plt.subplots(1, 1)
axs.plot(x, y)


x, y = runge_kutta(0.4, y_seta_ex, 2, [0, 5])

axs.plot(x, y)

x, y = runge_kutta(0.3, y_seta_ex, 2, [0, 5])

axs.plot(x, y)

x, y = runge_kutta(0.2, y_seta_ex, 2, [0, 5])

axs.plot(x, y)

x, y = runge_kutta(0.1, y_seta_ex, 2, [0, 5])

axs.plot(x, y)

plt.show()

# grafico(x, y)
