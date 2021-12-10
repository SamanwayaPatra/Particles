def solverRK4(self, Layout, idealz=0, idealvz=2, dt=0.01): 
        logger.info('Solver-RK4 Invoked')
        Nelts = len(Layout)
        MAXITER = 100
        logger.info(f'\t Layout is {[i.__str__() for i in Layout]}')

        # TODO: Simulate and find exit time rather than brue calculation
        t0 = 0
        idealts = np.empty(Nelts,dtype=np.float64)
        idealts[0] = t0 + (Layout[0].zf - idealz)/idealvz
        for i in range(1, Nelts):
            idealts[i] = idealts[i-1] + Layout[i].L/idealvz
        logger.info(f"Ideal Crossing times are {idealts}")


        tf = idealts[Nelts-1]
        tf = (Layout[-1].zi - Layout[0].zf)/idealvz
        logger.info(f'Calculated tf is {tf}')

        N = int((tf-t0)/dt) + 1
        logger.info(f'Number of Time Steps is : {N}')
        X0 = self.Xi
        t = np.zeros(N)
        X = np.zeros((N,6))
        Xcross = np.zeros((Nelts, 6))
        logger.info("Memory alloted for t, X, Xcross arrays")
        
        t[0] = t0
        X[0] = X0


        current_elt = 0
        tindex = 0
        xcross_index = 0

        while current_elt < Nelts:
            elt = Layout[current_elt]
            logger.info(f"Currently at Element: {current_elt}/{Nelts},\t: {elt}")
            observationtime = idealts[xcross_index]
            
            if elt.__str__() == 'Drift':
                (x, y, z, vx, vy, vz) = X[tindex]
                ti = t[tindex]
                delta_t = (elt.zf - z)/vz 
                tnext = t + delta_t 
                xnext = x + vx*delta_t
                ynext = y + vy*delta_t
                znext = z + vz*delta_t
                Xnext = np.array[xnext, ynext, znext, vx, vy, vz]
                tindex +=1 
                t[tindex] = tnext 
                X[i+1] = Xnext
                current_elt +=1

            else:
                insideElt = True 
                logger.info(f"Particle Inside {elt} ")
                while insideElt:
                    Xi = X[tindex]
                    ti = t[tindex]
                    logger.info(f"\tParticle's Config\t: {ti}, {Xi}")
                    tnext, Xnext = self.RK4Step(Xi, ti, dt, elt.field)
                    logger.info(f"\tParticle's Updated Config\t: {tnext}, {Xnext}")
                    tindex += 1
                    t[tindex] = tnext 
                    X[tindex+1] = Xnext
                    znext = Xnext[3]
                    
                    if tindex == observationtime:
                        self.Xcross[xcross_index] = Xnext
                        xcross_index +=1 
                     
                    if znext >= elt.zf:
                        insideElt = False 
                        current_elt +=1
                        break 

            if  xcross_index >= Nelts:
                break 

            if tindex >= MAXITER:
                break 