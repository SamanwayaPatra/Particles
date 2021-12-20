from numpy import float64
from particles.particle import Particle 
from particles.conversions import *
from particles.base_logger import logger
import matplotlib.pyplot as plt
plt.rcParams.update({'figure.autolayout': True})
class Particles():
    """ 
    Handle a beam of particle and its behaviour 

   Attributes
    ----------
    N - Number of Particles
    particles - An array of N particle objects

    Methods
    -------
    propagate particles - Solves the equation of motion for each particle

    """ 

    def __init__(self,Layout,observz,no_of_particles = 100): 
        self.No_of_Particles = no_of_particles
        self.t0 = 0
        self.Layout = Layout
        self.observz=observz
        self.Nelts = len(observz)
        #logger.info(f"Particles Class Instantiated with Layout {[i.__str__() for i in Layout]}")

    def readBeamParams(self, fname = "Beam_parameters.txt"):
        beam_data = np.loadtxt(fname)
        logger.info(f"Reading BEAM Data from File ({fname})")
        CHST, mass_number = beam_data[0], beam_data[1]
        Freq_MHz, no_of_particles = beam_data[2], int(beam_data[3])
        T_MeVpu, distribution_kind = beam_data[4], int(beam_data[5])
        a_x,  b_x,  e_x = beam_data[6],  beam_data[7],  beam_data[8]
        g_x = (1 + a_x*a_x)/b_x  
        a_y,  b_y,  e_y = beam_data[9],  beam_data[10], beam_data[11]       
        g_y = (1 + a_y*a_y)/b_y
        a_z,  b_z,  e_z = beam_data[12], beam_data[13], beam_data[14]  
        g_z = (1 + a_z*a_z)/b_z
        # TODO: Conversion of Parameters
        
        self.N = no_of_particles
        self.EMeVpu = T_MeVpu
        self.gamma = gamma_from_energy_MeVpu(T_MeVpu)
        self.beta = beta_from_gamma(self.gamma)
        self.mass_number = mass_number
        self.mass_kg = mass_number*amu_to_kg
        self.relativisticmass_kg = self.mass_kg*self.gamma
        self.momentum_Gevpc = momentum_from_energy_MeVpu_and_massnumber(T_MeVpu, mass_number)
        self.freq_Hz = Freq_MHz*10**6

        # TODO: What is this??
        self.wavelength = c_mps*10**3/self.freq_Hz
        
        self.TWISS = np.array(
                [   [a_x, b_x, g_x, e_x],
                    [a_y, b_y, g_y, e_y],
                    [a_z, b_z, g_z, e_z]
                ])
        logger.info(f"Sucessfully Read Beam Information from File")

    def isInsideEllipse(self, x, xdash, alpha, beta, gamma, epsilon):
        res = gamma*x*x + 2*alpha*x*xdash + beta*xdash*xdash - epsilon
        if res < 0 :
            return 1
        else :
            return 0
                     
    def twiss(self, TWIS=None):
        if TWIS == None : TWIS=self.TWISS
        logger.info(f"Calculating Twiss Parameters....")
        alpha_x, beta_x, gamma_x, epsilon_x = TWIS[0,0],  TWIS[0,1],   TWIS[0,2],   TWIS[0,3]
        alpha_y, beta_y, gamma_y, epsilon_y = TWIS[1,0],  TWIS[1,1],   TWIS[1,2],   TWIS[1,3]
        alpha_z, beta_z, gamma_z, epsilon_z = TWIS[2,0],  TWIS[2,1],   TWIS[2,2],   TWIS[2,3]
        logger.info(f"Converting Twiss Parameters to Co-ordinate Limits....")
        x_max = np.sqrt(beta_x*epsilon_x)
        y_max = np.sqrt(beta_y*epsilon_y)
        z_max = np.sqrt(beta_z*epsilon_z)
        xdash_max = np.sqrt(gamma_x*epsilon_x)
        ydash_max = np.sqrt(gamma_y*epsilon_y)
        zdash_max = np.sqrt(gamma_z*epsilon_z)
        
        self.MAXBEAMPARAMS = {"x_max":x_max, "y_max":y_max, "z_max":z_max, 
                           "xdash_max":xdash_max, "ydash_max":ydash_max, "zdash_max":zdash_max }
        self.TWISSPARAMS = {"alpha_x":alpha_x, "alpha_y":alpha_y, "alpha_z":alpha_z, 
                            "beta_x":beta_x, "beta_y":beta_y, "beta_z":beta_z, 
                            "gamma_x":gamma_x, "gamma_y":gamma_y, "gamma_z":gamma_z,
                            "epsilon_x":epsilon_x, "epsilon_y":epsilon_y, "epsilon_z":epsilon_z, 
                             }
        logger.info(f"TWISS Paramters\n {self.TWISSPARAMS}")
        logger.info(f"Coordinate Limiting Params\n {self.MAXBEAMPARAMS}")

    def generateUniformBeam(self):
        Xi = np.empty(self.N)
        Xdashi = np.empty(self.N)
        Yi = np.empty(self.N)
        Ydashi = np.empty(self.N)
        Zi = np.empty(self.N)
        Zdashi = np.empty(self.N)

        x_max, y_max, z_max = self.MAXBEAMPARAMS["x_max"], self.MAXBEAMPARAMS["y_max"], self.MAXBEAMPARAMS["z_max"]
        xdash_max, ydash_max, zdash_max = self.MAXBEAMPARAMS["xdash_max"], self.MAXBEAMPARAMS["ydash_max"], self.MAXBEAMPARAMS["zdash_max"]
        alpha_x, alpha_y, alpha_z = self.TWISSPARAMS["alpha_x"], self.TWISSPARAMS["alpha_y"], self.TWISSPARAMS["alpha_z"]
        beta_x, beta_y, beta_z = self.TWISSPARAMS["beta_x"], self.TWISSPARAMS["beta_y"], self.TWISSPARAMS["beta_z"]
        gamma_x, gamma_y, gamma_z = self.TWISSPARAMS["gamma_x"], self.TWISSPARAMS["gamma_y"], self.TWISSPARAMS["gamma_z"]
        epsilon_x, epsilon_y, epsilon_z = self.TWISSPARAMS["epsilon_x"], self.TWISSPARAMS["epsilon_y"], self.TWISSPARAMS["epsilon_z"]

        no_of_particles_generated = 0
        while no_of_particles_generated < self.N:
            x, y, z = np.random.uniform(-x_max, x_max), np.random.uniform(-y_max, y_max), np.random.uniform(-z_max, z_max)
            xdash, ydash, zdash = np.random.uniform(-xdash_max, xdash_max), np.random.uniform(-ydash_max, ydash_max), np.random.uniform(-zdash_max, zdash_max)
            is_in_x = self.isInsideEllipse(x, xdash, alpha_x, beta_x, gamma_x, epsilon_x)
            is_in_y = self.isInsideEllipse(y, ydash, alpha_y, beta_y, gamma_y, epsilon_y)
            is_in_z = self.isInsideEllipse(z, zdash, alpha_z, beta_z, gamma_z, epsilon_z)
            i = no_of_particles_generated
            if  is_in_x and is_in_y and is_in_z:
                Xi[i], Yi[i], Zi[i] = x, y, z 
                Xdashi[i], Ydashi[i], Zdashi[i] = xdash, ydash, zdash
                no_of_particles_generated = no_of_particles_generated + 1
        
        return (Xi, Xdashi, Yi, Ydashi, Zi, Zdashi )

    def convertBeamtoParticles(self, BEAM = None):
        if BEAM is None: BEAM = self.generateUniformBeam()
        logger.info("Creating Particle Distribution from BEAM data")
        (Xi, Xdashi, Yi, Ydashi, Zi, Zdashi ) = BEAM 
        self.particles = np.empty(self.N, dtype=type(Particle(self.Layout)))
        X1 = np.empty(self.N, dtype=np.float64); Y = np.empty(self.N, dtype=np.float64); Z = np.empty(self.N, dtype=np.float64); Vx = np.empty(self.N, dtype=np.float64); Vy = np.empty(self.N, dtype=np.float64); Vz = np.empty(self.N, dtype=np.float64)
        pz_Gevpc =  self.momentum_Gevpc
        PZ = pz_Gevpc*GeVpc_to_kgmps
        logger.info(f"Average Momentum of Particles is {PZ}")
        mass = self.relativisticmass_kg
        logger.info(f"Average Velocity of Beam {self.beta*c_mps}")
        logger.info("-"*10 + "Generating Particles from BEAM" + "-"*10)
        for i in range(self.N):
            pz = PZ + Zdashi[i]*PZ
            vz = pz/mass
            px = Xdashi[i]*pz
            vx = px/mass
            py = Ydashi[i]*pz
            vy = py/mass 

            x = Xi[i]
            y = Yi[i]
            z = Zi[i]
            InitialX = (x, y, z, vx, vy, vz)
            self.particles[i] = Particle(InitialX)
            # logger.info(f"Generated Particle with Initial Coodrinates {InitialX}")
            self.idealvz = PZ/mass
            self.idealz = 0
        logger.info("-"*10 + f"Generated {self.N} Particles from BEAM" + "-"*10)
        logger.info("Generated {self.N} Particles form BEAM")
        X1[i] = x; Y[i] = y; Z[i] = z; Vx[i] = vx; Vy[i] = vy; Vz[i] = vz 
        return (X1, Y, Z, Vx, Vy, Vz)
    
    def plotBeamPhase(self, BEAM = None):
        if BEAM is None: BEAM = self.generateUniformBeam()
        (Xi, Xdashi, Yi, Ydashi, Zi, Zdashi ) = BEAM 
        fig, axs = plt.subplots(3, figsize = (10,10))
        axs[0].scatter(Xi, Xdashi)
        axs[0].set_xlabel("x"); axs[0].set_ylabel("x'")
        axs[1].scatter(Yi, Ydashi)
        axs[1].set_xlabel("y"); axs[1].set_ylabel("y'")
        axs[2].scatter(Zi, Zdashi)
        axs[2].set_xlabel("z"); axs[2].set_ylabel("z'")
        # fig.show()

    def plotParticlePhase(self):
        (X, Vx, Y, Vy, Z, Vz) = self.convertBeamtoParticles()
        fig, axs = plt.subplots(3, figsize = (10,10))
        axs[0].scatter(X, Vx)
        axs[0].set_xlabel("x"); axs[0].set_ylabel("Vx")
        axs[1].scatter(Y, Vy)
        axs[1].set_xlabel("y"); axs[1].set_ylabel("Vy")
        axs[2].scatter(Z, Vz)
        axs[2].set_xlabel("z"); axs[2].set_ylabel("Vz")
        # fig.show()
        

    def createUniformDistribution(self,params=None): 
        if params is None: params = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 , 2 ]
        xmin, xmax, ymin, ymax, zmin, zmax, vxmin, vxmax, vymin, vymax, vzmin, vzmax = params 
        self.Xi = np.random.uniform(xmin, xmax, self.N)
        self.Yi = np.random.uniform(ymin, ymax, self.N)
        self.Zi = np.random.uniform(zmin, zmax, self.N)
        self.VXi = np.random.uniform(vxmin, vxmax, self.N)
        self.VYi = np.random.uniform(vymin, vymax, self.N)
        self.VZi = np.random.uniform(vzmin, vzmax, self.N)
        X = zip(self.Xi, self.Yi, self.Zi, self.VXi, self.VYi, self.VZi)
        self.particles = np.empty(self.N, dtype=type(Particle(self.Layout)))
        self.Xcross=np.zeros((self.N,self.Nelts,6),dtype=np.float64)
        for i,Xi in enumerate(X):
            self.particles[i] = Particle(Xi)
        
        self.idealz=(zmin+zmax)/2
        self.idealvz=(vzmax+vzmin)/2
    
    def idealts(self,obsz): #TODO change when implementing Electric Fields
        retts=np.empty([len(obsz)],dtype=float64)
        for i in range(len(obsz)):
            retts[i]=obsz[i]/self.idealvz
        return retts

    def propagateParticles(self): 
        logger.info("Simulating Paricle Trajectory")
        obsts=self.idealts(self.observz)
        self.X=np.empty([self.N,self.Nelts,6],dtype=float64)
        for i in range(self.N):
            print("Calculating Trajectory for Particle : {}/{}".format(i+1, self.N))
            logger.info("Calculating Trajectory for Particle : {}/{}".format(i+1, self.N))
            self.X[i]=self.particles[i].solverRK4(self.Layout,obsts)
                       
    def calcSig(self): 
        Nelts = self.Nelts
        Xideal = np.zeros( (Nelts, 6), np.float64 )
        deltaX = np.zeros( (Nelts, 6), np.float64 )
        X = np.zeros( (Nelts, 6), np.float64 )
        sigma = np.zeros( (Nelts,6,6), np.float64 )
        self.phase = np.zeros( (self.N, self.Nelts, 6), np.float64 )
        
        for i in range(Nelts):
            Xideal[i][2] = self.observz[i]
            Xideal[i][5] = self.idealvz
            
        for i in range(self.N):
            deltaX = np.subtract(self.X[i], Xideal)

            for j in range(Nelts):
                X[j][0] = deltaX[j][0]
                X[j][1] = deltaX[j][3]/self.idealvz
                X[j][2] = deltaX[j][1]
                X[j][3] = deltaX[j][4]/self.idealvz
                X[j][4] = deltaX[j][2]
                X[j][5] = deltaX[j][5]/self.idealvz 
            self.phase[i] = np.copy(X)
            for j in range(Nelts):
                for k in range(6):
                    for l in range(6):
                        sigma[j][k][l]+=(X[j][k]*X[j][l])/self.N
        return sigma

    def writePhase(self): 
        import xlwt
        wb=xlwt.Workbook()
        labels=["x","th(x)","y","th(y)","l","delta"]
        for idx in range(self.Nelts):
            sheet=wb.add_sheet("z="+str(self.observz[idx]))
            for i in range(6):
                sheet.write(0,i,labels[i])
            for i in range(self.N):
                for j in range(6):
                    sheet.write(i+1,j,self.phase[i][idx][j])
        wb.save(str(self.N)+"_"+str(self.Nelts)+".xls")

    def plotPhase(self,idx=0,xax=0,yax=1): 
        labels=["x","th(x)","y","th(y)","l","delta"]
        xaxis=self.phase[:,idx,xax]
        yaxis=self.phase[:,idx,yax]
        plt.plot(xaxis,yaxis,'r.')
        plt.xlabel(labels[xax])
        plt.ylabel(labels[yax])
        plt.title("Ideal Particle at z="+str(self.observz[idx]))
        plt.grid()
        plt.show()
    
    def plotPhases(self, idx=0):
        x = self.phase[:,idx,0]
        xdash = self.phase[:,idx,1]
        y = self.phase[:,idx,2]
        ydash = self.phase[:,idx,3]
        z = self.phase[:,idx,4]
        zdash = self.phase[:,idx,5]
        
        plt.figure(figsize = (10,7))
        
        plt.subplot(2,2,1)
        plt.title("Phase Space at z="+str(self.observz[idx]))
        plt.plot(x,xdash,'r.')
        plt.grid()
        plt.xlabel("x")
        plt.ylabel("x'")
        plt.subplot(2,2,2)
        plt.title("Phase Space at z="+str(self.observz[idx]))
        plt.plot(y,ydash,'r.')
        plt.grid()
        plt.xlabel("y")
        plt.ylabel("y'")
        plt.subplot(2,2,3)
        plt.title("Phase Space at z="+str(self.observz[idx]))
        plt.plot(z,zdash,'r.')
        plt.grid()
        plt.xlabel("z")
        plt.ylabel("z'")
