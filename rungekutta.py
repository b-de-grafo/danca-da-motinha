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
            y[k+1] = y[k] + (h/2) * (y_seta(x[k], y[k]) + y_seta(x[k] + h, y[k]+h*y_seta(x[k], y[k])))
        return x, y
    
    elif n > 1:
        y = [[y_zero[i]] for i in range(n)] # [[1000],[300]]
        print("y: " + str(y))
        x = list(np.arange(intervalo[0], intervalo[1], h))
        print("x: " + str(x))
        
        for k in range(len(x)-1):
            A = y_seta(x[k], [y[k] for y_k in y]) # ERRADO
            for i in range(n):
                print("i = " + str(i))
                y[i].append(None)
                print(f"y[i] = {y[i]}")
                
                print(f"y[i][k] = {y[i][k]}")
                # print(f"y_seta(x[{k}], y[{i}][{k}])[{i}] = {y_seta(x[k], y[i][k])[i]}")
                # dslkfds = y_seta(x[k] + h, y[i][k]+h*y_seta(x[k], y[i][k])[i])[i]
                # print(f"y_seta(x[{k}] + {h}, y[{i}][{k}]+{h}*y_seta(x[{k}], y[{i}][{k}])[{i}] = {dslkfds}")
                
                B = y[i][k] + h/2 * A
                C = x[k] + h/2
                D = y_seta(C, B)[i]
                print(f"A: {A}")
                print(f"B: {B}")
                print(f"C: {C}")
                print(f"D: {D}")
                
                # y[i][k+1] = y[i][k] + (h/2) * (y_seta(x[k], y[i][k])[i] + y_seta(x[k] + h, y[i][k]+h*y_seta(x[k], y[i][k])[i])[i])
                y[i][k+1] = y[i][k] + h * y_seta(x[k]+(h/2), y[i][k]+(h/2)*y_seta(x[k], y[i][k])[i])[i]
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

def lotka_volterra(x, y, a=100, b=0.37, c=100, d=0.05):
    print(f"y na lotka: {y}")
    return (a*y[0] - b*y[0]*y[1], -c*y[1] + d*y[0]*y[1])


fig, axs = plt.subplots(1, 1)

"""
x1, y1 = runge_kutta(0.1, y_seta_ex, 2, [0, 5])
axs.plot(x1, y1, color = "red")

x2, y2 = passo_multiplo(0.1, y_seta_ex, 2, [0, 5], runge_kutta)
axs.plot(x2, y2, color = "blue")
"""

x1, y1 = runge_kutta(h=0.002, y_seta=lotka_volterra, y_zero=[1000, 300], intervalo=[0,0.002*4], n=2)
print(x1)
print(y1)

plt.plot(x1, y1[0], color = "blue")
plt.plot(x1, y1[1], color = "green")

# plt.show()