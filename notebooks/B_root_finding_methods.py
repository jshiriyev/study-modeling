import numpy as np

import time

from scipy import optimize

class bisection():

    def __init__(self,func,lower,upper,tol=1e-5,args=None):

        self.lower = lower
        self.upper = upper

        self.tol = tol

        if args is None:
            args = []

        self._iterate(func,args)

    def _iterate(self,func,args):

        lower,upper = float(self.lower),float(self.upper)

        flower = func(lower,*args)
        fupper = func(upper,*args)

        self.calls = 2

        while (upper-lower)>self.tol:

            middle = (lower+upper)/2

            fmiddle = func(middle,*args)

            # print(f"{self.calls:2} {lower:8.6f},{flower:8.6f},{upper:8.6f},{fupper:8.6f}")

            self.calls += 1

            if fmiddle*flower<0:
                upper,fupper = middle,fmiddle
            elif fmiddle*fupper<0:
                lower,flower = middle,fmiddle
            else:
                self.info = "Either there is no root or even number of roots."
                middle,fmiddle = float('nan'),float('nan')
                break
        else:

            self.info = f"Found a root after {self.calls} iterations."

        self.value = middle

        self.minima = fmiddle

class newton():

    def __init__(self,func,der1,initial,tol=1e-5,args=None):

        self.initial = initial

        self.tol = tol

        self.calls = 0

        if args is None:
            args = []

        self._iterate(func,der1,args)

    def _iterate(self,func,der1,args):

        initial = self.initial

        d0 = func(initial,*args)
        d1 = der1(initial,*args)

        self.calls += 1

        while abs(d0/d1)>self.tol:

            print(f"{self.calls:2} {initial:8.6f},{d0:8.6f},{d1:8.6f},{d0/d1:8.6f}")

            initial -= d0/d1

            d0 = func(initial,*args)
            d1 = der1(initial,*args)

            self.calls += 1

        self.info = f"Converged after {self.calls} calls."

        self.value = initial

        self.minima = d0

class secant():

    def __init__(self,func,x0,x1,tol=1e-5,args=None):

        self.x0 = x0
        self.x1 = x1

        self.tol = tol

        self.calls = 0

        if args is None:
            args = []

        self._iterate(func,args)

    def _iterate(self,func,args):

        x0,x1 = self.x0,self.x1

        y0 = func(x0,*args)

        self.calls += 1

        # print(f"{self.calls:2} {x0:8.6f},{y0:8.6f}")

        y1 = func(x1,*args)

        self.calls += 1

        d1 = (y1-y0)/(x1-x0)

        while abs(y1/d1)>self.tol:

            # print(f"{self.calls:2} {x1:8.6f},{y1:8.6f},{y1/d1:8.6f}")

            x2 = x1-y1/d1

            y2 = func(x2,*args)

            self.calls += 1

            x0,y0 = x1,y1
            x1,y1 = x2,y2

            d1 = (y1-y0)/(x1-x0)

        self.value = x1

        self.minima = y1

if __name__ == "__main__":

    def f1(x):
        return x**2-4*np.sin(x)

    def f2(x):
        return 2*x-4*np.cos(x)

    sol1 = newton(f1,f2,3,tol=1e-8)

    def t1(x):
        return x**2-2*x+1

    def t2(x):
        return 2*x-2

    sol2 = newton(t1,t2,3,tol=1e-5)

    # sol3 = secant(f1,1,3,tol=1e-8)

    print(sol2.info)

    def obj(phi,epd,Re):
        xx = np.power(phi,-0.5)
        return xx+2*np.log10(epd/3.7+2.51/Re*xx)

    def obj_der(phi,epd,Re):
        xx = np.power(phi,-0.5)
        xy = np.power(phi,-1.5)
        return xy*(-1/2-2.51/np.log(10)/Re/(epd/3.7+2.51/Re*xx))

    epd = 1e-4

    Re = 100000

    phi0 = 8/Re
    # phi1 = 16/Re

    # tic = time.perf_counter()
    # phiS = root.secant(obj,phi0,phi1,tol=1e-7,args=(epd,Re)).value
    # toc = time.perf_counter()
    # print(f"home built secant method took \t {toc - tic:0.7f} seconds")

    tic = time.perf_counter()
    mysol = root.newton(obj,obj_der,phi0,tol=1e-7,args=(epd,Re))
    phiN = mysol.value
    toc = time.perf_counter()
    print(f"home built newton method took \t {toc - tic:0.7f} seconds")

    # tic = time.perf_counter()
    # phiM = optimize.newton(obj,phi0,args=(epd,Re))
    # toc = time.perf_counter()
    # print(f"scipy newton method took \t {toc - tic:0.7f} seconds\n")

    # print("home built secant method result is \t",phiS)
    # print("home built newton method result is \t",phiN)
    # print("scipy newton method result is \t\t",phiM)
    ##print(obj(0.01851386607747164,epd,Re))
    ##print(obj(0.01851386607747165,epd,Re))

    # print(mysol.info)



    def objective1(x):
        return x**3-x-2

    def objective2(x):
        return x**2-4*np.sin(x)

    def objective3(x):
        return x**2-4*x+4

    sol1 = root.bisection(objective1,1,2,tol=1e-5)

    print(sol1.info)
    print(sol1.value)
    print(sol1.minima,end="\n\n")

    sol2 = root.bisection(objective2,1,3,tol=1e-7)

    print(sol2.info)
    print(sol2.value)
    print(sol2.minima,end="\n\n")

    sol3 = root.bisection(objective3,1,3,tol=1e-7)

    print(sol3.info)
    print(sol3.value)
    print(sol3.minima)