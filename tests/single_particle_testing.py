import unittest
from particles.layout import trailLayout

class TestSingleParticle(unittest.TestCase):
    def test_solver(self):
        layout = trailLayout()
        print(layout)
