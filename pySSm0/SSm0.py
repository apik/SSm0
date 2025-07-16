import ctypes as ct
libSSm0 = ct.CDLL("libmmaSSm0.so")
libSSm0.SSm0_tldI.restype = ct.c_double
libSSm0.SSm0_tldS.restype = ct.c_double

def tldI(n, x, y):
    if x > 1 or x < 0 or y < -1 or y > 1 or n < -3 or n > 1 :
        return None
    else:
        return libSSm0.SSm0_tldI(n, ct.c_double(x), ct.c_double(y))

def tldS(n, x, y):
    if x > 1 or x < 0 or y < -1 or y > 1  or n < -3 or n > 1 :
        return None
    else:
        return libSSm0.SSm0_tldS(n, ct.c_double(x), ct.c_double(y))
