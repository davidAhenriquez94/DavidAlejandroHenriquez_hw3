import numpy as np

CO = np.zeros(9358)
NMHC  = np.zeros(9358)
C6H6 = np.zeros(9358)
NOx = np.zeros(9358)
NO2 = np.zeros(9358)

infile = open('AirQualityUCI.csv','r')
f = infile.readlines()
for i in range(1,9358):
    x = f[i].split(';')
    CO[i-1] = float(x[2].replace(',','.'))
    NMHC[i-1] = float(x[4].replace(',','.'))
    C6H6[i-1] = float(x[5].replace(',','.'))
    NOx[i-1] = float(x[7].replace(',','.'))
    NO2[i-1] = float(x[9].replace(',','.'))


N = 9358
a = 24
def periodo(n):
    return int(n*2*np.pi*(1.0/N)*a)
def exp(n):
    return np.exp(-2*np.pi*1j*(1.0/N)*n*a*np.arange(1,N+1))
    
def magnitud (z):
    return np.sqrt(np.real(z)**2 + np.imag(z)**2)
def fourier ():
    A = np.zeros((5,2))
    F = np.zeros((5,N+1))
    for n in range(N+1):
        F[0,n]=magnitud((1.0/a)*(np.sqrt((2*np.pi)))*np.dot(exp(n),CO))
        F[1,n]=magnitud((1.0/a)*(np.sqrt((2*np.pi)))*np.dot(exp(n),NMHC))
        F[2,n]=magnitud((1.0/a)*(np.sqrt((2*np.pi)))*np.dot(exp(n),C6H6))
        F[3,n]=magnitud((1.0/a)*(np.sqrt((2*np.pi)))*np.dot(exp(n),NOx))
        F[4,n]=magnitud((1.0/a)*(np.sqrt((2*np.pi)))*np.dot(exp(n),NO2))
    for i in range(5):
        A[i,0] = periodo(np.argmax(F[i,1::])+1)
        F[i,int(np.argmax(F[i,1::]))+1] = -100000
        A[i,1] = periodo(np.argmax(F[i,1::])+1)
    
    fileout  = open("periodos.dat", "w")

    for i in range(5):
        fileout.write("%f %f\n"%(A[i,0], A[i,1]))
    fileout.close()
    
fourier()

