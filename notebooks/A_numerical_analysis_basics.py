import numpy as np

def machine_epsilon(float_type):
    """Returns the number of significant digits in the mantissa
    for a binary floats.
    """

    Eps = float_type(1)
    One = float_type(1)
    Two = float_type(2)

    while True:
        
        if Eps/Two+One<=One:
            break
        
        Eps = Eps/Two

    t = (One-np.log(Eps)/np.log(Two))

    return int(t)

print(machine_epsilon(float))
print(machine_epsilon(np.float32))
print(machine_epsilon(np.float64))


# FLOATING POINT ARITHMETIC

# The objective of this exercise is to investigate the effect of round-off error
# on large numbers of interdependent computations.

sum1 = np.float32(0)    # single precision
sum2 = np.float64(0)    # double precision
sum3 = 0.               # double precision

x1 = np.float32(1e-5)
x2 = np.float64(1e-5)
x3 = 1e-5

for i in range(100000):
    sum1 += x1
    sum2 += x2
    sum3 += x3

print(sum1,abs(sum1-1)*100)
print(sum2,abs(sum2-1)*100)
print(sum3,abs(sum3-1)*100)


# CANCELLATION

def quadratic_equation(x,a=1,b=1,c=1):
    """Returns the value of quadratic equation."""    
    return a*x**2+b*x+c

def naive_solver(a,b,c):
    """Return the root values of x for quadratic equation,
    naive implementation."""

    D = np.sqrt(b**2-4*a*c)

    x1 = (-b+D)/(2*a)
    x2 = (-b-D)/(2*a)
        
    y1 = quadratic_equation(x1,a=a,b=b,c=c)
    y2 = quadratic_equation(x2,a=a,b=b,c=c)

    print(f"{x1=}, {y1=}")
    print(f"{x2=}, {y2=}")

def wise_solver(a,b,c):
    """Return the root values of x for quadratic equation,
    wise implementation."""

    D = np.sqrt(b**2-4*a*c)

    if np.sign(b*D) == 0:
        x1 = 0 if a!=0 else np.nan
        x2 = 0 if a!=0 else np.nan
    elif np.sign(b*D)>0:
        x1 = (2*c)/(-b-D)
        x2 = (-b-D)/(2*a)
    elif np.sign(b*D)<0:
        x1 = (-b+D)/(2*a)
        x2 = (2*c)/(-b+D)

    y1 = quadratic_equation(x1,a=a,b=b,c=c)
    y2 = quadratic_equation(x2,a=a,b=b,c=c)

    print(f"{x1=}, {y1=}")
    print(f"{x2=}, {y2=}")

if __name__ == "__main__":
    
    a = 1.
    b = 200.
    c = -1.5e-14

    naive_solver(a,b,c)
    wise_solver(a,b,c)

