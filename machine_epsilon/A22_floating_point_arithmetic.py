import numpy as np

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
