import numpy as np
import matplotlib.pyplot as plt
import math


def diferencas_finitas_exp(c_quadrado, delta_t, delta_x, intervalo, u_x_zero, cond_contorno, cond_contorno_derivada, passos_t):
    gama_quadrado = c_quadrado * (delta_t**2) / (delta_x**2)
    u = []

    # t0
    subu = []
    for x in np.arange(intervalo[0], intervalo[1]+delta_x, delta_x):
        if x == intervalo[0]:
            subu.append(cond_contorno[0](x, 0))
        elif x == intervalo[1]:
            subu.append(cond_contorno[1](x, 0))
        else:
            subu.append(u_x_zero(x, 0))
    u.append(subu)

    atualiza_grafico(np.arange(intervalo[0], intervalo[1]+delta_x, delta_x), subu)

    # t1
    subu = []
    i = 0
    for x in np.arange(intervalo[0], intervalo[1]+delta_x, delta_x):
        subu.append(delta_t*cond_contorno_derivada(x, delta_t) + u[0][i])
        i = i+1
    u.append(subu)

    atualiza_grafico(np.arange(intervalo[0], intervalo[1]+delta_x, delta_x), subu)

    # t > 1
    k = 1
    for t in np.arange(delta_t * 2, delta_t * passos_t, delta_t):
        subu = []
        i = 0
        for x in np.arange(intervalo[0], intervalo[1]+delta_x, delta_x):
            if x == intervalo[0]:
                u_i_k_mais_um = cond_contorno[0](0, t)
            elif x == intervalo[1]:
                u_i_k_mais_um = cond_contorno[1](0, t)
            else:
                u_i_k_mais_um = gama_quadrado * (u[k][i-1] - 2*u[k][i] + u[k][i+1]) + 2 * u[k][i] - u[k-1][i]
            subu.append(u_i_k_mais_um)
            i = i + 1
        u.append(subu)
        k = k + 1

        atualiza_grafico(np.arange(intervalo[0], intervalo[1] + delta_x, delta_x), subu)
    return u


def atualiza_grafico(x, u):
    ax1.clear()
    ax1.set_ylim(bottom=-1.5, top=1.5)
    line,  = ax1.plot(x, u)
    line.set_xdata(x)
    line.set_ydata(u)
    plt.draw()
    plt.pause(0.001)


def x_vezes_x_menos_1(x, t):
    return x*(x-1)


def sen_pi_vezes_t(x, t):
    return math.sin(math.pi*t)


def zero(x, t):
    return 0


def dois(x, t):
    return 2


fig = plt.figure()
ax1 = fig.add_subplot(1, 1, 1)
print(diferencas_finitas_exp(4, 0.01, 0.2, [0, 1], x_vezes_x_menos_1, [sen_pi_vezes_t, zero], dois, 10000))
