from particles.layout import trailLayout
from particles.particles import Particles
from particles.base_logger import logger
import numpy as np
Layout = trailLayout()
observz=np.array([0,0.1,0.2,0.7,1.2,1.3,1.4])
p = Particles(Layout,observz)
p.readBeamParams("data/Beam_parameters.txt")
p.twiss()
p.convertBeamtoParticles()
p.propagateParticles()
p.calcSig()
p.plotPhase(1)
p.writePhase()