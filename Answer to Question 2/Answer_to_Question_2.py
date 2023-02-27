# 8.2

import numpy as np
import matplotlib.pyplot as plt

# define initial value
t = 0.01  # s
T = 30  # s

T_ls = np.arange(0, T, t)
E_ls = np.zeros(T_ls.shape, float)
S_ls = np.zeros(T_ls.shape, float)
ES_ls = np.zeros(T_ls.shape, float)
P_ls = np.zeros(T_ls.shape, float)

E_ls[0] = 1  # µM
S_ls[0] = 10  # µM
ES_ls[0] = 0  # µM
P_ls[0] = 0  # µM

k1 = 100/60  # µM/s
k2 = 600/60  # µM/s
k3 = 150/60  # µM/s

# Define the equations
def func_E(E, S, ES):
    fE = -k1*E*S + k2*ES + k3*ES
    return fE

def func_S(E, S, ES):
    fS = -k1*E*S + k2*ES
    return fS

def func_ES(E, S, ES):
    fES = k1*E*S - k2*ES - k3*ES
    return fES

def func_P(ES):
    fP = k3*ES
    return fP

# The fourth-order RungeKutta method
for i in range(T_ls.size-1):
    E1 = func_E(E_ls[i], S_ls[i], ES_ls[i])
    S1 = func_S(E_ls[i], S_ls[i], ES_ls[i])
    ES1 = func_ES(E_ls[i], S_ls[i], ES_ls[i])
    P1 = func_P(ES_ls[i])

    E2 = func_E(E_ls[i]+E1*t/2, S_ls[i]+S1*t/2, ES_ls[i]+ES1*t/2)
    S2 = func_S(E_ls[i]+E1*t/2, S_ls[i]+S1*t/2, ES_ls[i]+ES1*t/2)
    ES2 = func_ES(E_ls[i]+E1*t/2, S_ls[i]+S1*t/2, ES_ls[i]+ES1*t/2)
    P2 = func_P(ES_ls[i]+ES1*t/2)

    E3 = func_E(E_ls[i]+E2*t/2, S_ls[i]+S2*t/2, ES_ls[i]+ES2*t/2)
    S3 = func_S(E_ls[i]+E2*t/2, S_ls[i]+S2*t/2, ES_ls[i]+ES2*t/2)
    ES3 = func_ES(E_ls[i]+E2*t/2, S_ls[i]+S2*t/2, ES_ls[i]+ES2*t/2)
    P3 = func_P(ES_ls[i]+ES2*t/2)

    E4 = func_E(E_ls[i]+E3*t/2, S_ls[i]+S3*t/2, ES_ls[i]+ES3*t/2)
    S4 = func_S(E_ls[i]+E3*t/2, S_ls[i]+S3*t/2, ES_ls[i]+ES3*t/2)
    ES4 = func_ES(E_ls[i]+E3*t/2, S_ls[i]+S3*t/2, ES_ls[i]+ES3*t/2)
    P4 = func_P(ES_ls[i]+ES3*t/2)

    E_ls[i+1] = E_ls[i] + t*(E1+2*E2+2*E3+E4)/6
    S_ls[i+1] = S_ls[i] + t*(S1+2*S2+2*S3+S4)/6
    ES_ls[i+1] = ES_ls[i] + t*(ES1+2*ES2+2*ES3+ES4)/6
    P_ls[i+1] = P_ls[i] + t*(P1+2*P2+2*P3+P4)/6

# Print the answer
print('8.2')
print('The final concentration of enzyme E is %.5f µM.' %E_ls[-1])
print('The final concentration of substrate S is %.5f µM.' %S_ls[-1])
print('The final concentration of intermediate species ES is %.5f µM.' %ES_ls[-1])
print('The final concentration of product P is %.5f µM.' %P_ls[-1])

# Draw the picture
plt.figure()
plt.plot(T_ls, E_ls, label='E(enzyme)')
plt.plot(T_ls, S_ls, label='S(substrate)')
plt.plot(T_ls, ES_ls, label='ES(compound)')
plt.plot(T_ls, P_ls, label='P(product)')
plt.xlabel('Time(s)')
plt.ylabel('Concentration(µM)')
plt.legend()
plt.show()


# 8.3

V_ls = k3*ES_ls
Vm = max(V_ls)*60
index = np.argmax(V_ls)
Sm = S_ls[index]

#Print the answer
print('8.3')
print('When the substrate concentration is %.5f µM, the reaction rate reaches its maximum value.'%Sm)
print ('The maximum speed Vm is %.5f µM/min.'%Vm)

#Draw the picture
plt.figure()
plt.plot(S_ls, V_ls)
plt.xlabel('Concentration of Substrate)(µM)')
plt.ylabel('Velocity(s)')
plt.show()

