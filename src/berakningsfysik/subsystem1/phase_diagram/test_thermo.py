import unittest
import math
from thermo import calculate_delta_s_mix, calculate_delta_g_mix

class TestThermodynamics(unittest.TestCase):
    def test_ideal_mixing(self):
            x = [0.5, 0.5] #blandning av två komponenter
            R = 8.314462618 #molar gas constant (J/(mol*K))
            expected_s = -R * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
            
            calculated_s = calculate_delta_s_mix(x, moles=1.0, R=8.314462618)
            self.assertAlmostEqual(calculated_s, expected_s, places=4)
    
    def test_delta_g_mix_temperature(self):
        x = [0.5, 0.5]
        T = 300
        S = calculate_delta_s_mix(x, moles=1.0, R=8.314462618)
        expected_g = -T * S
        
        calculated_g = - calculate_delta_g_mix(x, T, moles=1.0, R=8.314462618)
        self.assertAlmostEqual(calculated_g, expected_g, places=4)

    def test_pure_substance(self):
        x = [1.0, 0.0]
        self.assertEqual(calculate_delta_s_mix(x, moles=1.0, R=8.314462618), 0.0)
        #Equal to 0 because there is no mixing, only one component
    
    def test_not_one_mole_fractions(self):
        with self.assertRaises(ValueError):
            calculate_delta_s_mix([0.5, 0.2], moles=1.0, R=8.314462618) #Sum != 1

    def test_negative_mole_fractions(self):
        with self.assertRaises(ValueError):
            calculate_delta_s_mix([-0.1, 1.1], moles=1.0, R=8.314462618) #There is a negative number

if __name__ == '__main__':
    unittest.main()

#S = Entropy
#G = Gibbs free energy
#H = Enthalpy
#Delta S = R * sum(x_i * ln(x_i)) for i in components
#Delta G = Delta H - T * Delta S
#Delta H = 0 for ideal mixing, so Delta G = -T * Delta S