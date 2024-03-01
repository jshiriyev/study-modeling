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

    import numpy as np

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