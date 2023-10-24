import numpy,time

def h2(dist=0.75):
    if dist == 0.7:     return numpy.load("matrix/H2_0.70.npy")
    if dist == 0.75:    return numpy.load("matrix/H2_0.75.npy")
    if dist == 1.5:     return numpy.load("matrix/H2_1.5.npy")
    raise ValueError("Unsupported bond length. Use `legacymolecule` instead?")

def lih(dist=1.5):
    if dist == 1.5:     return numpy.load("matrix/LiH_1.5.npy")
    if dist == 3.0:     return numpy.load("matrix/LiH_3.0.npy")
    raise ValueError("Unsupported bond length. Use `legacymolecule` instead?")