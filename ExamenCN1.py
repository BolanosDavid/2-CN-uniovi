"""
Ejercicios Primer Control Computacion Numerica
"""
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import numpy.polynomial.polynomial as pol
# Funciones usadas en cada practica
"""Practica4"""
def horner(x0, p):
    q = np.zeros_like(p)
   
    for i in range(len(p)-1, -1, -1):
       if i == len(p)-1:
           q[i] = p[i]
       else:
           q[i] = p[i] + q[i+1]*x0
   
    cociente = q[1:]
    resto = q[0]
    return cociente, resto
def HornerV(x, p):
    y = np.zeros_like(x)
    for j in range(len(x)):
        q = np.zeros_like(p)
        for i in range(len(p)-1,-1,-1):
            if(i == len(p)-1):
                q[i] = p[i]
            else:
                q[i] = p[i] + q[i+1]*x[j]
        y[j] = q[0]
    return y
def dersuc(p,x0):
    n = len(p)
    der = np.zeros_like(p)
    restos = np.zeros_like(p)
    factorial = 1.
    for i in range(n):
        q, resto = horner(x0,p)
        restos[i] = resto
        der[i]=resto*factorial
        p=q
        factorial *=i+1
    return der,restos
def divisores(m):
    if m == 0:
        return np.array([0])
    div_list = []
    for i in range(1, m + 1):
        if m % i == 0:
            div_list.append(i)
            div_list.append(-i)
    
    return np.array(div_list)
def Raices(p):
    p_copy = p.copy()
    raices = []
    
    while len(p_copy) > 1: 
        m = int(abs(p_copy[0]))
        if m == 0:
            raices.append(0)
            p_copy = p_copy[:-1] 
            continue
        posibles_raices = divisores(m)
        raiz_encontrada = False
        for r in posibles_raices:
            q, resto = horner(r, p_copy)
            if resto == 0:
                raices.append(r)
                p_copy = q
                raiz_encontrada = True
                break
        
        if not raiz_encontrada:
            break
    
    return np.array(raices)
      

"""Practica5"""
def imprimirFuncion(f, r, a, b, n=100):
    x = np.linspace(a, b, n)
    y = f(x) 
    ox = np.zeros_like(x) 
    plt.plot(x, y, label="f(x)")
    plt.plot(x, ox, 'k', label="Eje X")
    if len(r) > 0:
        plt.plot(r, np.zeros_like(r), 'ro', label="Raíces")
    plt.title('Gráfica de la función con raíces')
    plt.legend()
    plt.show()              
def bolzano(f,a,b):
    return f(a)*f(b) < 0
def busquedaIncremental(f,a,b,n=100):
    x = np.linspace(a,b,n+1)
    intervalos = np.zeros((n,2))
    y = f(x)
    contador = 0
    
    for i in range(n):
        if y[i]*y[i+1] < 0:
            intervalos[contador,:] = x[i:i+2]
            contador += 1
            
    intervalos = intervalos[:contador,:]    
    return intervalos
def biseccion(f,a,b,tol=1.e-6,n=100):
    if f(a)*f(b) > 0.0: 
        print('La función tiene el mismo signo en los extremos')
    
    for i in range(n):
        m = 0.5*(a + b)
                
        if ( bolzano(f,a,m)): 
            b = m;
        elif bolzano(f,m,b):  
            a = m
        else:
            break
        
        if b - a < tol:
            break
    
    return m,i+1
def newton(f,df,x0,tol=1.e-6,n=100):
    for i in range(n):
        x1 = x0 - f(x0)/df(x0)
        if abs(x1 - x0) < tol: 
                            
            break
        x0 = x1
    return x1,i+1
def raices_bisec(f,a,b):
    intervalos = busquedaIncremental(f, a, b)
    raices = []      
    iteraciones = [] 
    for xa, xb in intervalos:
        r, i = biseccion(f, xa, xb)
        raices.append(r)
        iteraciones.append(i)

    return np.array(raices), np.array(iteraciones)
def raices_newton(f,df,a,b):
    intervalos = busquedaIncremental(f, a, b)
    raices = []
    iteraciones = []
    for xa,xb in intervalos:
        r , i = newton(f, df, xa)
        raices.append(r)
        iteraciones.append(i)
    return raices,iteraciones    
def secante(f,x0,x1,tol = 1.e-6,maxiter=100):
    for it in range(maxiter):
        x2 = x1-f(x1)*((x1-x0)/(f(x1)-f(x0)))
        if(abs(x2 - x1) < tol):
            break
        x0 = x1
        x1 = x2
    return x2,it+1
def punto_fijo(g, x0, tol = 1E-6, max_iter = 100):
    x = [x0]
    for i in range(max_iter):
        x_k = g(x[-1])
        x.append(x_k)
        if np.abs(x[-1] - x[-2]) < tol:
            return x[-1], x, i
    return x[-1], x, max_iter
"""Ejercicios extra"""
def cambios_signo(p):
    cont = 0
    for i in range(0, len(p) - 1):
        ceros = 0  
        while i + ceros + 1 < len(p) and p[i + ceros + 1] == 0:
            ceros += 1
        if i + ceros + 1 < len(p) and p[i] * p[i + ceros + 1] < 0:
            cont += 1
        i += ceros
    return cont
def derivada(p):
    p_Derivado = np.zeros(len(p)-1)
    counter = 0
    for i in range(len(p)-1, 1, -1):
        p_Derivado[counter]=p[counter]*i
        counter += 1
        
    return p_Derivado

"""Practica4"""
#%% Ejercicio1 Horner:ruffini
   
p0 = np.array([1., 2, 1])
x0 = 1.
rp0 =  pol.polyval(x0,p0) 

p1 = np.array([1., -1, 2, -3,  5, -2])
x1 = 1.
rp1 =  pol.polyval(x1,p1) 

p2 = np.array([1., -1, -1, 1, -1, 0, -1, 1])
x2 = -1.
rp2 =  pol.polyval(x2,p2) 


c, r = horner(x0, p0)
print("Cociente = ", c)
print("Resto = ", r)
print('Con polyval = ', rp0)

c1, r1 = horner(x1, p1)
print("Cociente = ", c1)
print("Resto = ", r1)
print('Con polyval = ', rp1)

c2, r2 = horner(x2, p2)
print("Cociente = ", c2)
print("Resto = ", r2)  
print('Con polyval = ', rp2)    
 
#%%Ejercicio2 HornerV(x,p)
x = np.linspace(-1, 1)
p1 = np.array([1., -1, 2, -3, 5, -2])
r = HornerV(x,p1)
plt.figure()
plt.plot(x, r)
plt.plot(x, x*0, 'k', label="Eje X")
plt.title('Gráfico del Polinomio1 usando el Método de Horner')
plt.show()
p2 = np.array([1., -1, -1, 1, -1, 0, -1, 1])
r = HornerV(x,p2)
plt.figure()
plt.plot(x, r)
plt.plot(x, x*0, 'k', label="Eje X")
plt.title('Gráfico del Polinomio2 usando el Método de Horner')
plt.show()
#%%Ejercicio3 dersuc(x0,p)
np.set_printoptions(suppress = True)
p1 = np.array([1., -1, 2, -3,  5, -2])
x0 = 1.
derivadas,restos = dersuc(p1, x0)
print("Restos de dividir P1 una y otra vez por (x-",int(x0),")")
print(restos)
print("Derivadas sucesivas de P en x0 = ", x0)
print(derivadas)

p2 = np.array([1., -1, -1, 1, -1, 0, -1, 1])
x0 = -1.
derivadas,restos = dersuc(p2, x0)
print("Restos de dividir P1 una y otra vez por (x-",int( x0),")")
print(restos)
print("Derivadas sucesivas de P en x0 = ", x0)
print(derivadas)
#%% Ejercicio4  raices()
print("Divisores de 6")
print(divisores(6))
print("Divisores de 18")
print(divisores(18))
print("Divisores de 20")
print(divisores(20))
p0 = np.array([-1.,0,1])
p1 = np.array([8., -6, -3, 1])
p2 = np.array([15., -2, -16, 2, 1])
p3 = np.array([60.,53, -13, -5, 1])   
p4 = np.array([490., 343, -206, -56, 4, 1])


print("Raices de p0 = ",Raices(p0))
print("Raices de p1 = ",Raices(p1))
print("Raices de p2 = ",Raices(p2))
print("Raices de p3 = ",Raices(p3))
print("Raices de p4 = ",Raices(p4))
#%% Ejercicio5 multiples raices
p1 = np.array([8., -22, 17, 1, -5, 1])
p2 = np.array([-135., 378, -369, 140, -9, -6, 1])
p3 = np.array([96., 320, 366, 135, -30, -24, 0, 1]) 
p4 = np.array([280., 156, -350, -59, 148, -26, -6, 1])
print("Raices de p1 = ",Raices(p1))
print("Raices de p2 = ",Raices(p2))
print("Raices de p3 = ",Raices(p3))
print("Raices de p4 = ",Raices(p4))

"""Practica5"""    
#%% Ejercicio1 busquedaIncremental
"""Ejemplo1"""
f1 = lambda x: x**5-3*x**2 + 1.6
a,b,n = -1,1.5,25
print(busquedaIncremental(f1,a,b,n))
"""Ejemplo2"""
f1 = lambda x: (x+2)*np.cos(2*x)
a,b,n = 0,10,100
print(busquedaIncremental(f1,a,b,n))

#%% Ejercicio2 biseccion
""" Ejemplo 1"""
f1 = lambda x: x**5 - 3 * x**2 + 1.6 
a_original = -1; b_original = 1.5; n = 25
intervalos = busquedaIncremental(f1,a_original,b_original,n)
r = np.zeros(len(intervalos))
i = np.zeros(len(intervalos))
for x in range(r.size):
    a,b = intervalos[x]
    r[x],i[x] = biseccion(f1,a,b)
imprimirFuncion(f1, r, a_original, b_original)
print("Raices obtenidas mediante biseccion:\n",r,"\nIteraciones utilizadas: ",i)
""" Ejemplo 2"""
f1 = lambda x: ((x**3+1)/(x**2+1))*np.cos(x) - 0.2
a = -3; b = 3; n = 100
intervalos = busquedaIncremental(f1,a,b,n)
r = np.zeros(len(intervalos))
i = np.zeros(len(intervalos))
for x in range(r.size):
    a,b = intervalos[x]
    r[x],i[x] = biseccion(f1,a,b)
imprimirFuncion(f1, r, a_original, b_original)
print("Raices obtenidas mediante biseccion:\n",r,"\nIteraciones utilizadas: ",i)
#%% Ejercicio3 newton

x = sym.Symbol('x', real=True)
"""Ejemplo1"""
f_sim   = x**5-3*x**2+1.6
f = sym.lambdify(x, f_sim)
df_sim  = sym.diff(f_sim,x)
df = sym.lambdify(x, df_sim) 
a = -1.; b = 1.5; n = 25
intervalos = busquedaIncremental(f,a,b,n)
r = np.zeros(len(intervalos))
i = np.zeros(len(intervalos))
for x in range(r.size):
    a,b = intervalos[x]
    r[x],i[x] = newton(f,df,a)
print("Raices obtenidas mediante newton:\n",r,"\nIteraciones utilizadas: ",i)
imprimirFuncion(f, r, a, b)
"""Ejemplo 2"""
x = sym.Symbol('x', real=True)
f_sim   = ((x**3+1)/(x**2+1))*sym.cos(x) - 0.2
f = sym.lambdify(x, f_sim)
df_sim  = sym.diff(f_sim,x)
df = sym.lambdify(x, df_sim) 
a = -3.; b = 3; n = 25
intervalos = busquedaIncremental(f,a,b,n)
r = np.zeros(len(intervalos))
i = np.zeros(len(intervalos))
for x in range(r.size):
    a,b = intervalos[x]
    r[x],i[x] = newton(f,df,a)
print("Raices obtenidas mediante newton:\n",r,"\nIteraciones utilizadas: ",i)
imprimirFuncion(f, r, a, b)
#%% Ejercicio4 raices_bisec
x = sym.Symbol('x', real=True)
f_sim = x**5 - 3*x**2 + 1.6 
f = sym.lambdify(x, f_sim, 'numpy')  # Convertir la función simbólica a función de NumPy
a, b = -2, 2
raices, iteraciones = raices_bisec(f, a, b)
print("Raíces encontradas:", raices)
print("Iteraciones usadas por cada raíz:", iteraciones)
imprimirFuncion(f, raices, a, b)
#%%Ejercicio5 raices_newton
"""Ejemplo1"""
x = sym.Symbol('x', real=True)
f_sim = x**4 + 2*x**3 - 7*x**2 + 3
df_sim = sym.diff(f_sim)
df = sym.lambdify(x,df_sim)
f = sym.lambdify(x, f_sim)
a, b = -4, 2
raices, iteraciones = raices_newton(f,df, a, b)
imprimirFuncion(f, raices, a, b)
"""Ejemplo2"""
x = sym.Symbol('x', real=True)
f_sim = x**6 - 0.1*x**5 - 17*x**4 + x**3 + 73*x**2 - 4 * x - 68
df_sim = sym.diff(f_sim)
df = sym.lambdify(x,df_sim)
f = sym.lambdify(x, f_sim)
a, b = -4, 4
raices, iteraciones = raices_newton(f,df, a, b)
imprimirFuncion(f, raices, a, b)
#%%Ejercicio6 secante
np.set_printoptions(precision = 5)
f = lambda x: x**5 -3*x**2 + 1.6
a,b,n = -1,1.5,25
intervalos = busquedaIncremental(f, a, b,n)
r = np.zeros(3)
iteraciones = np.zeros(3)
for x in range(3):
    a,b = intervalos[x]
    r[x],iteraciones[x] = secante(f, a, b)
print(r,iteraciones)
""""Ejercicios extras"""
#%% Ejercicio cambios de signo
p0 = np.array([32., -32, -14, 17, -3])
p1 = np.array([32., -32, 0, 17, -3])
T4 = np.array([8, 0, -8, 0, 1])
T6 = np.array([32., 0, -48, 0, 18, 0, -1])

print("P0")
print(cambios_signo(p0), " raíces reales positivas como máximo")
print("P1")
print(cambios_signo(p1), " raíces reales positivas como máximo")
print("T4")
print(cambios_signo(T4), " raíces reales positivas como máximo")
print("T6")
print(cambios_signo(T6), " raíces reales positivas como máximo")