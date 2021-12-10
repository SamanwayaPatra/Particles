import numpy as np 
from particles.base_logger import logger

class Particle():
    """ 
    Handle a single particle and its behaviour 

    Attributes
    ----------
    Q - Charge of Particle 
    m - Mass of Particle 
    Xi - Initial position-velocity vector 
    traj - Contains the entire trajectory of the particle as an object containg tuple of (t,X) t is timeslices and X is position velocity vector at each timestep. Created when SolverRK4 is run.



    Methods
    -------
    dXdt(X,t, field) 
        Takes a vector X containing position and velocity (x,y,z,vx,vy,vz) and returns its derivative which is (vx,vy,vz,ax,ay,az)

    solverInbuilt(t0, tf, X0)
        Takes a time range (t0 and tf) with an initial vector X0 and solves the equation of motion (dXdt) using scipy's odeint method
    
    RK4Step(t0,X0,dt)
        Takes an initial vector at X0 at t0 and returns the new vector at time t0 + dt under the DE dXdt using RK4 algoritm.
        Called only by solver RK4. 
    
    solverRK4(t0,tf,dt)
        Solves the equation of motion of particle for a time range (t0, tf) using the initial condition Xi for the equation of motion dXdt using RK4 Algorithm for the timestep dt
        Returns traj object which is a tuple of (t,X) t -> array containing timesteps X -> array containing position veloctiy vector for each timestep also saves this into to traj in class attribute

    """


    def __init__(self, Xi=None, Q=1, m=1): 
        if Xi is None: Xi =[0, 0, 0, 0, 0, 0]
        self.Xi = Xi
        self.Q = Q 
        self.m = m
        logger.info(f"Particle Created with Initial Condition X = {Xi}")
       
    # ------------------ TO BE COMPLETED --------------- #
    # TODO: Complete later
    def lorentzForce():
        pass

    def dXdt(self, X, t, field):
        x,y,z,vx,vy,vz = X
        Q = self.Q 
        m = self.m 
        Bx, By, Bz = field(x, y, z)
        # TODO: Allow for Electric Forces
        dvxdt = np.array(Q/m*(vy*Bz - vz*By), dtype=np.float64)
        dvydt = np.array(Q/m*(vz*Bx - vx*Bz), dtype=np.float64)
        dvzdt = np.array(Q/m*(vx*By - vy*Bx), dtype=np.float64)
        dxdt = vx 
        dydt = vy 
        dzdt = vz 
        dXdt = np.array((dxdt, dydt, dzdt, dvxdt, dvydt, dvzdt))
        return dXdt
        
    def solverInbuilt(self, ti = 0, tf = 10, X0 = None):
        if X0 is None: X0 = [0, 0, 0, 0, 1, 2]
        from scipy.integrate import odeint
        t = np.linspace(ti,tf,100)
        sol = odeint(self.dXdt, X0, t, )
        return t,sol 

    def RK4Step(self, t0, X0, dt, field):
        Bx, By, Bz = field(X0[0], X0[1], X0[2])
        ODE = self.dXdt
        K1 = dt*np.array(ODE(X0, t0, field))
        K2 = dt*np.array(ODE(X0+0.5*K1,t0 + dt*0.5, field))
        K3 = dt*np.array(ODE(X0+0.5*K2,t0 + dt*0.5, field))
        K4 = dt*np.array(ODE(X0+K3, t0 + dt, field))
        X = X0 + (K1 + 2*K2 + 2*K3 + K4)/6
        X = np.array(X)
        t = t0 + dt
        logger.info(f"Change in X is: {[K1, K2, K3, K4]}")
        logger.info(f"Particle's Updates Config is \t({t}, {X}) under field of ({Bx}, {By}, {Bz})") 
        return (t,X)

    def solverRK4(self, Layout, idealz=0, idealvz=2, ti=0, dt=1e-4): 
        Nelts = len(Layout)
        logger.info(f"RK-4 Solver invoked for Particle in Layout: {[i.__str__() for i in Layout]}")

        # TODO: Simulate and find exit time rather than brue calculation
        # TODO: Make this Global
        idealts = np.empty(Nelts,dtype=np.float64)
        idealts[0] = ti + (Layout[0].zf - idealz)/idealvz
        for i in range(1, Nelts):
            idealts[i] = idealts[i] + Layout[i-1].L/idealvz
        tf = 10


        N = int((tf-ti)/dt) + 1
        tf = idealts[Nelts-1]
        X0 = self.Xi
        t = np.zeros(N)
        X = np.zeros((N,6))
        self.Xcross=np.zeros((len(Layout),6),dtype=np.float64)
        t[0] = ti
        X[0] = X0
        tnew = ti
        Eltidx = i = j = 0
        logger.info(f"No of time steps is {N}, Memory sucesfully allocated for arrays: X, t")
        logger.info(f"Particle's Initial Config is \t({ti}, {X0})")
        
        # TODO: Redo this Parts....
        while tnew < tf :
            ti = t[i]
            Xi = X[i]
            tnew, Xnew = self.RK4Step(ti,Xi,dt,Layout[Eltidx].field)
            logger.info(f"Particle's Updates Config is \t({tnew}, {Xnew}) under field of {Layout[Eltidx].__str__()}")
            
            if tnew >= idealts[j]:
                self.Xcross[j] = Xnew
                j+=1
                if j == Nelts:
                    break

            if Xnew[2] >= Layout[Eltidx].zf :
                Eltidx+=1
                if Eltidx == Nelts:
                    while j<Nelts:
                        tnew=idealts[j]-tnew
                        Xnew[0]+=tnew*Xnew[3]
                        Xnew[1]+=tnew*Xnew[4]
                        Xnew[2]+=tnew*Xnew[5]
                        tnew=idealts[j]
                        self.Xcross[j]=Xnew
                        j+=1
                    break
                elif Layout[Eltidx].__str__() == "Drift":
                    X[i+1] = Xnew
                    t[i+1] = tnew
                    i+=1
                    tnew=(Layout[Eltidx].zf-Xnew[2])/Xnew[5]
                    Xnew[0]+=tnew*Xnew[3]
                    Xnew[1]+=tnew*Xnew[4]
                    Xnew[2]+=tnew*Xnew[5]
                    tnew+=t[i]
                    while idealts[j]<tnew:
                        self.Xcross[j][0]=X[i-1][0]+(idealts[j]-t[i-1])*X[i-1][3]
                        self.Xcross[j][1]=X[i-1][1]+(idealts[j]-t[i-1])*X[i-1][4]
                        self.Xcross[j][2]=X[i-1][2]+(idealts[j]-t[i-1])*X[i-1][5]
                        self.Xcross[j][3]=X[i-1][3]
                        self.Xcross[j][4]=X[i-1][4]
                        self.Xcross[j][5]=X[i-1][5]
                        j+=1
                        if j==Nelts:
                            return
                    Eltidx+=1
                    if Eltidx==Nelts:
                        while j<Nelts:
                            tnew=idealts[j]-tnew
                            Xnew[0]+=tnew*Xnew[3]
                            Xnew[1]+=tnew*Xnew[4]
                            Xnew[2]+=tnew*Xnew[5]
                            tnew=idealts[j]
                            self.Xcross[j]=Xnew
                            j+=1
                        break
            X[i+1] = Xnew
            t[i+1] = tnew
            i+=1
        #self.traj = Trajactory(t,X)
        return (t,X)
