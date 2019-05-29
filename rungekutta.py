import numpy as np
import matplotlib.pyplot as plt


def runge_kutta(h, y_seta, y_zero, intervalo):
    y = [y_zero]
    x = list(np.arange(intervalo[0], intervalo[1], h))
    for k in range(len(x)-1):
        y.append(None)
        y[k+1] = y[k] + (h/2) * (y_seta(x[k], y[k]) + y_seta(x[k] + h, y[k]+h*y_seta(x[k], y[k])))
    return x, y


def runge_kutta(h, y_seta, y_zero, intervalo, n=1):
    if n == 1:
        y = [y_zero]
        x = list(np.arange(intervalo[0], intervalo[1], h))
        for k in range(len(x)-1):
            y.append(None)
            y[k+1] = y[k] + h * y_seta(x[k] + h/2, y[k] + h/2 * y_seta(x[k],y[k]))
        return x, y
    
    elif n > 1:
        y = [[y_zero[i] for i in range(n)]]  # [[1000, 300]]
        x = list(np.arange(intervalo[0], intervalo[1], h))
        b = []
        for k in range(len(x)-1):
            a = y_seta(x[k], y[k])
            y.append([])
            b.append([])
            for i in range(n):
                y[k+1].append(None)
                b[k].append(y[k][i] + h/2 * a[i])
            c = x[k] + h/2
            d = y_seta(c, b[k])
            for i in range(n):
                y[k + 1][i] = y[k][i] + h * d[i]

        return x, y
        

def passo_multiplo(h, y_seta, y_zero, intervalo, inicializacao):
    y = [y_zero]
    x = list(np.arange(intervalo[0], intervalo[1], h))
    
    y.append(inicializacao(h, y_seta, y_zero, [intervalo[0], intervalo[0]+2*h])[1][1])    
    
    for k in range(1, len(x)-1):
        y.append(None)
        y[k+1] = y[k] + (h/2) * (3 * y_seta(x[k], y[k]) - y_seta(x[k-1], y[k-1]))
    
    return x, y


def passo_multiplo(h, y_seta, y_zero, intervalo, inicializacao, n=1):
    if n == 1:
        y = [y_zero]
        x = list(np.arange(intervalo[0], intervalo[1], h))

        y.append(inicializacao(h, y_seta, y_zero, [intervalo[0], intervalo[0] + 2 * h])[1][1])

        for k in range(1, len(x) - 1):
            y.append(None)
            y[k + 1] = y[k] + (h / 2) * (3 * y_seta(x[k], y[k]) - y_seta(x[k - 1], y[k - 1]))

        return x, y

    elif n > 1:
        y = [[y_zero[i] for i in range(n)]]  # [[1000, 300]]
        x = list(np.arange(intervalo[0], intervalo[1], h))

        y.append(inicializacao(h, y_seta, y_zero, [intervalo[0], intervalo[0] + 2 * h], n=n)[1][1])
        for k in range(1,len(x) - 1):
            a = y_seta(x[k], y[k])
            b = y_seta(x[k-1], y[k-1])
            y.append([])
            for i in range(n):
                y[k + 1].append(None)
                y[k + 1][i] = y[k][i] + h/2 * (3*a[i]-b[i])
        return x, y


def y_seta_ex(x, y):
    return -y+x-2


def lotka_volterra(x, y, a=100, b=0.37, c=100, d=0.05):
    return (a*y[0] - b*y[0]*y[1], -c*y[1] + d*y[0]*y[1])


def eq_motoca(x, y, k_1=21000, k_2=1500, m_1=40, m_2=400, b=1500, u=0):
    return (y[1], (k_2 * (y[2] - y[0]) + b * (y[3] - y[1]) + k_1 * (u - y[0])) / m_1, y[3], (-k_2 * (y[2] - y[0]) - b * (y[3] - y[2])) / m_2)


def transpoe_matriz(matriz):
    resp = [[matriz[j][i] for j in range(len(matriz))] for i in range(len(matriz[0]))]
    return resp




"""
x1, y1 = runge_kutta(0.1, y_seta_ex, 2, [0, 5])
axs.plot(x1, y1, color = "red")

x2, y2 = passo_multiplo(0.1, y_seta_ex, 2, [0, 5], runge_kutta)
axs.plot(x2, y2, color = "blue")
"""

"""
fig, axs = plt.subplots(2, 1)
x, y_rk = runge_kutta(h=0.002, y_seta=lotka_volterra, y_zero=[1000, 300], intervalo=[0, 0.002*100], n=2)
print(y_rk)
x, y_pm = passo_multiplo(h=0.002, y_seta=lotka_volterra, y_zero=[1000, 300], intervalo=[0, 0.002*100], inicializacao=runge_kutta, n=2)
print(y_pm)

y_rk = transpoe_matriz(y_rk)
print(y_rk)

y_pm = transpoe_matriz(y_pm)
print(y_pm)

axs[0].plot(x, y_rk[0], color="blue")
axs[0].plot(x, y_rk[1], color="blue")
axs[0].plot(x, y_pm[0], color="green")
axs[0].plot(x, y_pm[1], color="green")
axs[0].set_ylabel("populção")
axs[0].set_xlabel("tempo")

axs[1].plot(y_rk[0], y_rk[1], color="blue")
axs[1].plot(y_pm[0], y_pm[1], color="green")
axs[1].set_ylabel("predador")
axs[1].set_xlabel("presa")
plt.show()
"""

fig, axs = plt.subplots(1, 1)
x, y_rk = runge_kutta(h=0.01, y_seta=eq_motoca, y_zero=[0, 1, 0, 1], intervalo=[0, 10], n=4)
print(x)
print(y_rk)
x, y_pm = passo_multiplo(h=0.01, y_seta=eq_motoca, y_zero=[0, 1, 0, 1], intervalo=[0, 10], inicializacao=runge_kutta, n=4)
print(y_pm)

y_rk = transpoe_matriz(y_rk)
print(y_rk)

y_pm = transpoe_matriz(y_pm)
print(y_pm)

axs.plot(x, y_rk[0], color="blue")
axs.plot(x, y_rk[1], color="green")
axs.plot(x, y_rk[2], color="red")
axs.plot(x, y_rk[3], color="yellow")
axs.plot(x, y_pm[0], color="blue")
axs.plot(x, y_pm[1], color="green")
axs.plot(x, y_pm[2], color="red")
axs.plot(x, y_pm[3], color="yellow")
axs.set_xlabel("tempo")


plt.show()
