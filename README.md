![](https://raw.githubusercontent.com/rafafrdz/CC_P1/master/img/1.png)

![](https://raw.githubusercontent.com/rafafrdz/CC_P1/master/img/2.png)

```python
# -*- coding: utf-8 -*-
import time
from numpy import *
from scipy.sparse.linalg import cgs,splu
from scipy.sparse import lil_matrix,identity
from scipy.linalg import lu_factor,lu_solve,cho_factor,cho_solve
from matplotlib.pyplot import *


def f0(x):
    y=5.0*exp(-(x-0.5)**2)
    return y

#Resuelve el problema de contorno
# -nu*u'' +u=f
# u(x0)=ua
# u(xf)=ub
#usando matrices vacias y el resolvedor directo
#ver http://docs.scipy.org/doc/scipy/reference/sparse.html
#ver http://docs.scipy.org/doc/scipy/reference/sparse.linalg.html#module-scipy.sparse.linalg

####################################################
###EJERCICIO 1
def Ejercicio1(x0,xf,N,alfa,ua,fuente):
    t1=time.time()
    N=int(N)
    x0=float(x0)
    xf=float(xf)
    dx=(xf-x0)/float(N)
    dx2=dx*dx
    x=linspace(x0,xf,N+1)
    M = lil_matrix((N+1,N+1), dtype='float64')
    Id=identity(N+1,dtype='float64',format='csc')

 #    for i in range(0,N+1):
 #        M[i,i]=2.0
	# if (i<N):
	#    M[i,i+1]=-1.0
	# if (i>0):
	#    M[i,i-1]=-1.0
    M.setdiag(2.0*ones(N+1),0)
    M.setdiag(-1.0*ones(N+1),1)
    M.setdiag(-1.0*ones(N+1),-1)
    b=fuente(x)
    M[0,0]=0.0
    M[0,1]=0.0
    M[1,0]=0.0
    M=M.tocsc()
    A=Id+alfa/(dx2)*M
    A[N,N] = 1 + alfa/dx2
    A[N-1,N] = - alfa/dx2
    #imposicion de condiciones de contorno
    #modificacion del vector b
    b[0]=ua
    b[1]+=ua*alfa/dx2


    #hace la descomposición LU completa de una matriz Sparse
    LU=splu(A)
    #resolvemos el sistema lineal usando la descomposición LU
    usol=LU.solve(b)

    tf=time.time()
    print "Tiempo de ejecucion:",tf-t1
    plot(x,usol,'b')
    show()
#Ejercicio1(0,1,100,2,1,f0)

####################################################
###EJERCICIO 2
def u0(x):
    return 0*x+1
def Ejercicio2(x0,xf,N,T,Nt,alfa,ua,ub,fuente,cond0):
        t1=time.time()
        N=int(N)
        x0=float(x0)
        xf=float(xf)
        dx=(xf-x0)/float(N)
        dx2=dx*dx
        dt = 0.05
        x=linspace(x0,xf,N+1)
        u0 = cond0(x)
        M = lil_matrix((N+1,N+1), dtype='float64')
        Id=identity(N+1,dtype='float64',format='csc')

        M.setdiag(2.0*ones(N+1),0)
        M.setdiag(-1.0*ones(N+1),1)
        M.setdiag(-1.0*ones(N+1),-1)

        M[0,0]=0.0
        M[0,1]=0.0
        M[N,N]=0.0
        M[N,N-1]=0.0

        M[N-1,N]=0.0
        M[1,0]=0.0
        M=M.tocsc()
        A=Id+(alfa*dt)/(dx2)*M

        usol = u0
        LU=splu(A)
        for i in range(Nt):
            b=dt*fuente(x) + usol
            b[0]=ua
            b[N]=ub
            usol=LU.solve(b)
            plot(x,usol,'b')


        tf=time.time()
        print "Tiempo de ejecucion:",tf-t1
        show()
#Ejercicio2(0,1,100,2,20,2,0,0,f0,u0)


####################################################
###EJERCICIO3
def Ejercicio3(x0,xf,N,T,Nt,alfa,ua,fuente,cond0):
    t1=time.time()
    N=int(N)
    x0=float(x0)
    xf=float(xf)
    dx=(xf-x0)/float(N)
    dx2=dx*dx
    dt = T/float(Nt)
    x=linspace(x0,xf,N+1)
    M = lil_matrix((N+1,N+1), dtype='float64')
    Id=identity(N+1,dtype='float64',format='csc')

    M.setdiag(2.0*ones(N+1),0)
    M.setdiag(-1.0*ones(N+1),1)
    M.setdiag(-1.0*ones(N+1),-1)
    b=fuente(x)
    M[0,0]=0.0
    M[0,1] = 0.0
    M=M.tocsc()
    A=Id+(alfa*dt)/(dx2)*M
    A[N,N] = 1 + alfa*dt/dx2

    usol = cond0(x)
    LU=splu(A)
    for i in range(Nt+1):
        b=dt*fuente(x)+usol
        b[0]=ua
        usol=LU.solve(b)
        plot(x,usol,'b')


    tf=time.time()
    print "Tiempo de ejecucion:",tf-t1
    plot(x,usol,'b')
    show()
#Ejercicio3(0,1,100,2,40,2,1,f0,u0)
```

