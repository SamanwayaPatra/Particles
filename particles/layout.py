from numpy import float64
from particles.magnets import * 

def trailLayout():
    z = 0
    #q1 = Quadrapole(0, 0.2, 0.05, z)
    q1=Drift(0.2,z)
    z = q1.zf 
    d1 = Drift(1,z)
    z = d1.zf
    q2=Drift(0.2,z)
    #q2 = Quadrapole(0, 0.2, 0.05, z)
    z = q2.zf

    Lay=Layo([q1,d1,q2])
    return Lay

class Layo():
    def __init__(self,Lay):
        self.Lay=Lay
        self.nelts=len(Lay)
        self.bounds=np.zeros([self.nelts,2],dtype=float64)
        for i in range(self.nelts):
            self.bounds[i]=Lay[i].getBounds()
    def field(self,x,y,z):
        Bx=By=Bz=0
        for i in range(self.nelts):
            if(z>=self.bounds[i,0] and z<self.bounds[i,1]):
                B=self.Lay[i].field(x,y,z)
                Bx+=B[0]
                By+=B[1]
                Bz+=B[2]
        return np.array([Bx,By,Bz])
