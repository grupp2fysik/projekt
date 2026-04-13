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

    def test_different_mole_count(self):
        x = [0.5, 0.5]
        moles = 2.0
        R = 8.314462618
        expected_s = -R * moles * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        
        calculated_s = calculate_delta_s_mix(x, moles=moles, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=4)

    def test_more_than_one_component(self):
        x = [0.25, 0.25, 0.25, 0.25] #blandning av fyra komponenter
        R = 8.314462618
        expected_s = -R * 4 * (0.25 * math.log(0.25)) #Each component contributes equally to the entropy
        
        calculated_s = calculate_delta_s_mix(x, moles=1.0, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=4)
    
    def test_division_of_three(self):
        x = [0.333333333333333, 0.333333333333333, 0.333333333333333] #blandning av tre komponenter
        R = 8.314462618
        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, moles=1.0, R=R) #Sum != 1

    def test_pure_substance(self):
        x = [1.0, 0.0]
        R = 8.314462618
        self.assertEqual(calculate_delta_s_mix(x, moles=1.0, R=R), 0.0)
        #Equal to 0 because there is no mixing, only one component
    
    def test_not_one_mole_fractions(self):
        x = [0.5, 0.2]
        R = 8.314462618
        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, moles=1.0, R=R) #Sum != 1

    def test_negative_mole_fractions(self):
        x = [-0.1, 1.1]
        R = 8.314462618
        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, moles=1.0, R=R) #There is a negative number

    def test_invalid_mole_fractions(self):
        x = [0.5, 0.5, 0.5] #Sum > 1
        R = 8.314462618
        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, moles=1.0, R=R) #Sum != 1
    
    def test_negative_temperature(self):
        x = [0.5, 0.5]
        T = -50 #Negative temperature is not valid
        with self.assertRaises(ValueError):
            calculate_delta_g_mix(x, T, moles=1.0, R=8.314462618)

    def test_zero_components(self):
        x = [] #No components
        R = 8.314462618
        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, moles=1.0, R=R) #No components to mix
    
    def test_almost_one(self):
        x = [0.999999999999999, 1e-15] #Sum is very close to 1
        R = 8.314462618
        expected_s = -R * (0.999999999999999 * math.log(0.999999999999999) + 1e-15 * math.log(1e-15))
        
        calculated_s = calculate_delta_s_mix(x, moles=1.0, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=4)
    
    def test_negative_moles(self):
        x = [0.5, 0.5]
        moles = -1.0 #Negative moles is not valid
        R = 8.314462618
        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, moles=moles, R=R)
        
    def test_valid_gas_constant(self):
        x = [0.5, 0.5]
        R = 8.314462618 #Valid gas constant
        calculated_s = calculate_delta_s_mix(x, moles=1.0, R=R)
        self.assertIsNotNone(calculated_s)

    def test_zero_temperature(self):
        x = [0.5, 0.5]
        T = 0 #At absolute zero, the Gibbs free energy change should be zero
        R = 8.314462618
        expected_g = 0.0
        
        calculated_g = - calculate_delta_g_mix(x, T, moles=1.0, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=4)
    
    def test_zero_moles(self):
        x = [0.5, 0.5]
        moles = 0 #No moles means no mixing, so entropy change should be zero
        R = 8.314462618
        expected_s = 0.0
        
        calculated_s = calculate_delta_s_mix(x, moles=moles, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=4)
    
    def test_very_high_temperature(self):
        x = [0.5, 0.5]
        T = 1e6
        #At very high temperatures, the Gibbs free energy change should be dominated by the entropy term
        R = 8.314462618
        expected_g = -T * calculate_delta_s_mix(x, moles=1.0, R=R)
        
        calculated_g = - calculate_delta_g_mix(x, T, moles=1.0, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=4)

    def test_very_high_moles(self):
        x = [0.5, 0.5]
        moles = 1e6
        R = 8.314462618
        expected_s = -R * moles * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        
        calculated_s = calculate_delta_s_mix(x, moles=moles, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=4)
    
    def test_very_low_temperature(self):
        x = [0.5, 0.5]
        T = 1e-6
        R = 8.314462618
        expected_g = -T * calculate_delta_s_mix(x, moles=1.0, R=R)
        
        calculated_g = - calculate_delta_g_mix(x, T, moles=1.0, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=4)

    def test_very_low_moles(self):
        x = [0.5, 0.5]
        moles = 1e-6
        R = 8.314462618
        expected_s = -R * moles * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        
        calculated_s = calculate_delta_s_mix(x, moles=moles, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=4)
    
    def test_invalid_type(self):
        x = "not a list" #Invalid type for mole fractions
        R = 8.314462618
        with self.assertRaises(TypeError):
            calculate_delta_s_mix(x, moles=1.0, R=R)

if __name__ == '__main__':
    unittest.main()

#S = Entropy
#G = Gibbs free energy
#H = Enthalpy
#Delta S = R * sum(x_i * ln(x_i)) for i in components
#Delta G = Delta H - T * Delta S
#Delta H = 0 for ideal mixing, so Delta G = -T * Delta S