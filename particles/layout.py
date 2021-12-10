from particles.magnets import * 

def trailLayout():
    z = 0
    q1 = Quadrapole(20, 0.2, 0.05, z)
    z = q1.zf 
    d1 = Drift(1,z)
    z = d1.zf 
    q2 = Quadrapole(20, 0.2, 0.05, z)
    z = q2.zf

    Layout=[q1,d1,q2]
    return Layout