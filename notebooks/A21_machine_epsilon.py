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

if __name__ == "__main__":

    print(machine_epsilon(float))
    print(machine_epsilon(np.float32))
    print(machine_epsilon(np.float64))
