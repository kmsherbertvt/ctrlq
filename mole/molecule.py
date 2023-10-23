import numpy,time

def h2(dist=0.75):
    assert dist==0.75 or dist==1.5
    return numpy.load("matrix/h207.npy") if dist==0.75 else numpy.load("matrix/h215.npy")

def lih(dist=1.5):
    assert dist==1.5 or dist==3.0
    return numpy.load("matrix/lih15.npy") if dist == 1.5 else numpy.load("matrix/lih30.npy")