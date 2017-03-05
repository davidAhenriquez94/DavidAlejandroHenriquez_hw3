import numpy as np
import matplotlib.pyplot as plt
import os

filein = open('condiciones.dat','r')
datos =  filein.read().splitlines()
y_in = [int(i) for i in datos[4].split()]
v_in = [int(i) for i in datos[5].split()]

def leapfrog(N,T,m,a,t,h):
    m = int(t/h)
    Y = np.zeros((m,N+2))
    V = np.zeros((m,N+2))
    v_half_in = np.zeros(N+2)
    for i in range(m):
        for j in range(N+2):
            if i == 0:
                Y[i,1:N+1] = y_in
                V[i,1:N+1] = v_in
                v_half_in[1] = V[i,1] + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i,1]-Y[i,2])
                v_half_in[N] = V[i,N] + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i,N]-Y[i,N-1])
                if 1<j<N:
                    v_half_in[j] = V[i,j] + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i,j]-Y[i,j-1]-Y[i,j+1])
            
            else :
                if 0<j<N:
                    if i == 1:
                        Y[i,j] = Y[i-1,j] + h*v_half_in[j]
                        Y[i,j+1] = Y[i-1,j] + h*v_half_in[j+1]
                        V[i,j] = v_half_in[j] + (1/2.0)*h*(-T/m*a)*(2.0*Y[i,j]-Y[i,j-1]-Y[i,j+1])
                    else :
                        v_half_1 = V[i-1,j] + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i-1,j]-Y[i-1,j-1]-Y[i-1,j+1])
                        Y[i,j] = Y[i-1,j] + h*v_half_1
                        v_half_2 = V[i-1,j+1] + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i-1,j+1]-Y[i-1,j]-Y[i-1,j+2])
                        Y[i,j+1] = Y[i-1,j+1] + h*v_half_2
                        V[i,j] = v_half_1 + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i,j]-Y[i,j-1]-Y[i,j+1])
                elif j == N:
                    if i == 1:
                        Y[i,j] = Y[i-1,j] + h*v_half_in[j]
                        Y[i,j+1] = Y[i-1,j+1] + h*v_half_in[j+1]
                        V[i,j] = v_half_in[j] + (1/2.0)*h*(-T/m*a)*(2.0*Y[i,j]-Y[i,j-1]-Y[i,j+1])
                    else :
                        v_half = V[i-1,j] + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i-1,j]-Y[i-1,j-1]-Y[i-1,j+1])
                        Y[i,j] = Y[i-1,j] + h*v_half
                        V[i,j] = v_half + (1.0/2)*h*(-T/(m*a))*(2.0*Y[i,j]-Y[i,j-1]-Y[i,j+1])
                    
                    
    x = np.arange(N+2)*a
    for i in range(m):
        fig = plt.figure()
        ax = plt.axes()
        plt.plot(x,Y[i,::],'o-',label='N = '+str(N))
        plt.ylim(-np.amax(Y)*1.1,np.amax(Y)*1.1)
        ax.set_title("Simulacion Osciladores acoplados")
        ax.set_xlabel("tiempo")
        ax.set_ylabel("posicion")
        ax.legend()
        plt.savefig(str(i)+'.png', format = 'png')
        plt.close()
    os.system('convert -delay 50 -loop 0 $(ls *.png | sort -n) movimiento.gif')
    os.system('rm -f *.png')
    
    for i in range(m):
        os.system('rm -f ' + str(i)+'.gif')

    E = []
    for i in range(m):
        K = []
        U = []
        for j in range(1,N+1):
            K.append((1/2.0)*m*(V[i,j])**2)
            U.append((1/2.0)*(T/a)*(Y[i,j-1]-Y[i,j])**2)
        E.append(sum(K)+sum(U))
        
    fig = plt.figure()
    plt.plot((E-E[0])/E[0])
    plt.savefig('energia.pdf', format = 'pdf')
    plt.close()

leapfrog(int(datos[0]),float(datos[1]),float(datos[2]),float(datos[3]),200.0,0.5)
