from SSm0 import tldI, tldS
# This script tests the functions tldI and tldS from the SSm0 module.
# It prints the results for quarks and gluons with specified parameters.

epord = 1                       # epsilon order
x = 0.1                         # beta
y = 0.2                         # cos(theta)

# quarks
print("\\tilde{{I}}(x={},y={}, ep^{}) = {}".format(x,y,epord, tldI(epord, x ,y)))
# gluons
print("\\tilde{{S}}(x={},y={}, ep^{}) = {}".format(x,y,epord, tldS(epord, x ,y)))

