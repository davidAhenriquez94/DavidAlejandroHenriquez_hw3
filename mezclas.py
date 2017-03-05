import numpy as np
import matplotlib.pyplot as plt

infile = open('AirQualityUCI.csv','r')
f = infile.readlines()

N = 9357
A = np.zeros((N,5))
I = []
for i in range(N):
   A[i,0] = f[i+1].split(';')[2].replace(',','.')
   A[i,1] = f[i+1].split(';')[4].replace(',','.')
   A[i,2] = f[i+1].split(';')[5].replace(',','.')
   A[i,3] = f[i+1].split(';')[7].replace(',','.')
   A[i,4] = f[i+1].split(';')[9].replace(',','.')
   for j in range(5):
      if A[i,j] == -200:
         I.append(i)
         break

   A[i,1] = A[i,1]*(10.0**(-3))
   A[i,2] = A[i,2]*(10.0**(-3))
   A[i,3] = A[i,3]*1.64
   A[i,4] = A[i,4]*(10.0**(-3))

B = np.delete(A,I,axis=0)
N = B.shape[0] 

def cov_matrix(x):
   n = x.shape
   C = np.zeros((n[1],n[1]))
   for i in range(n[1]):
      for j in range(n[1]):
         C[i,j] = np.dot((x[::,i] - np.ones(n[0])*np.mean(x[::,i])),(x[::,j]-np.ones(n[0])*np.mean(x[::,j])))*(1.0/(N))
   return C

def comp_principales(X):
   y1,y2 = np.linalg.eig(cov_matrix(X))
   max_1 = np.argmax(y1)
   y_copy = np.copy(y1)
   y_copy[max_1] = -10000000
   max_2 = np.argmax(y_copy)
   return (y2[max_1,::],y2[max_2,::]),(max_1,max_2)

def grafica_PCA(X):
   x = []
   y = []
   for i in range(N):
      k1,k2 = comp_principales(X)
      x.append(np.dot(X[i,::],k1[0])/(np.sqrt(np.dot(X[i,::],X[i,::])))*(np.sqrt(np.dot(k1[0],k1[0]))))
      y.append(np.dot(X[i,::],k1[1])/(np.sqrt(np.dot(X[i,::],X[i,::])))*(np.sqrt(np.dot(k1[1],k1[1]))))
   
   fig = plt.figure()
   plt.scatter(x,y)
   plt.savefig('pca.pdf',format = 'pdf')
   plt.close()
   
def datos_PCA(X):
   fileout = open('pca.dat','w')
   k1,k2 = comp_principales(X)
   x = k1[0]
   y = k1[1]
   for i in range(len(x)):
      fileout.write("%f %f\n" %(x[i],y[i]) )
   fileout.close()

def varianzas_principales(X):
   M = cov_matrix(X)
   k1,k2 = comp_principales(X)
   i = int(k2[0])
   j = int(k2[1])
   
   sigm_1 = (M[i,i]/np.matrix.trace(M))*100
   sigm_2 = (M[j,j]/np.matrix.trace(M))*100
   
   fileout = open('varianza.dat','w')
   fileout.write("%f %f\n" %(sigm_1,sigm_2))

grafica_PCA(B)
datos_PCA(B)
varianzas_principales(B)

