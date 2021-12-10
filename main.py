from particles.layout import trailLayout
from particles.particles import Particles
from particles.base_logger import logger

Layout = trailLayout()
p = Particles(Layout)
p.readBeamParams("data/Beam_parameters.txt")
p.twiss()
p.convertBeamtoParticles()
p.propagateParticles()