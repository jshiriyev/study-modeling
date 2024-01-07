import numpy as np

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


