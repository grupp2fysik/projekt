import unittest
import math
from thermo import calculate_delta_s_mix, calculate_delta_g_mix, calculate_delta_g_mix_with_h

class TestDeltaSMix(unittest.TestCase):
    def test_ideal_mixing(self):
        x = [0.5, 0.5] #mixing of two components in equal proportions
        R = 8.314462618 #molar gas constant (J/(mol*K))
        n = 1.0 #number of moles

        expected_s = -R * (0.5 * math.log(0.5) + 0.5 * math.log(0.5)) #Expected entropy change for ideal mixing of two components in equal proportions
        calculated_s = calculate_delta_s_mix(x, n=n, R=R) #Assert that the calculated entropy change is approximately equal to the expected value, with a very high precision
        self.assertAlmostEqual(calculated_s, expected_s, places=100) #The places=100 argument allows for a very high precision in the comparison, which is important for testing the correctness of the function with very small differences in the expected and calculated values.

    def test_delta_s_bigger_than_zero(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0

        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertGreater(calculated_s, 0) #Assert that the calculated entropy change is greater than zero, which is expected for a mixing process that increases disorder.

    def test_different_mole_count(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 2.0

        expected_s = -R * n * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

    def test_more_than_one_component(self):
        x = [0.25, 0.25, 0.25, 0.25]
        R = 8.314462618
        n = 1.0

        expected_s = -R * 4 * (0.25 * math.log(0.25)) #Each component contributes equally to the entropy
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)
    
    def test_division_of_three(self):
        x = [0.333333333333333, 0.333333333333333, 0.333333333333333]
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #Sum != 1

    def test_pure_substance(self):
        x = [1.0, 0.0]
        R = 8.314462618
        n = 1.0

        self.assertEqual(calculate_delta_s_mix(x, n=n, R=R), 0.0)
        #Equal to 0 because there is no mixing, only one component
    
    def test_mole_fractions_not_equal_to_one(self):
        x = [0.5, 0.2] #Sum != 1
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #Sum != 1

    def test_negative_mole_fractions(self):
        x = [-0.1, 1.1] #There is a negative mole fraction, which is not valid
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)
    
    def test_invalid_mole_fractions(self):
        x = [0.5, 0.5, 0.5] #Sum > 1
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #Sum != 1

    def test_no_mole_fractions(self):
        x = [] #No components
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #No components to mix
    
    def test_almost_one(self):
        x = [0.9999999999999999, 1e-16] #Sum = 1, but one component is almost negligible
        R = 8.314462618
        n = 1.0

        expected_s = -R * (0.9999999999999999 * math.log(0.9999999999999999) + 1e-16 * math.log(1e-16))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)
 
    def test_negative_n(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = -1.0 #Negative n is not valid

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)

    def test_zero_n(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 0 #No n means no mixing, so entropy change should be zero

        expected_s = 0.0
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

    def test_very_high_n(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1e6 #Very high n means a very large entropy change

        expected_s = -R * n * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

    def test_very_low_n(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1e-6 #Very low n means a very small entropy change

        expected_s = -R * n * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

class TestDeltaGMix(unittest.TestCase):
    def test_delta_g_mix_temperature(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300 #Tempature in Kelvin

        expected_g = -T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_pure_substance(self):
        x = [1.0, 0.0]
        R = 8.314462618
        n = 1.0
        T = 300

        expected_g = 0.0 #No mixing, so Gibbs free energy change should be zero
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_negative_temperature(self):
        x = [0.5, 0.5]
        R = 8.314462618
        T = -50 #Negative temperature is not valid
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_g_mix(x, T, n=n, R=R)

    def test_zero_temperature(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 0 #At absolute zero, the Gibbs free energy change should be zero

        expected_g = 0.0
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_very_low_temperature(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e-6 #Very low temperature means a very small Gibbs free energy change

        expected_g = -T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_very_high_temperature(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e6 #Very high temperature means a very large Gibbs free energy change

        expected_g = -T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_temperature_relationship(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T1 = 300
        T2 = 600

        g1 = calculate_delta_g_mix(x, T1, n=n, R=R)
        g2 = calculate_delta_g_mix(x, T2, n=n, R=R)
        self.assertAlmostEqual(g2, 2 * g1, places=100)

class TestDeltaGMixWithH(unittest.TestCase):
    def test_delta_g_mix_with_delta_h(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 1000 #Enthalpy change for mixing

        expected_g = delta_h - T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_delta_g_mix_with_delta_h_pure_substance(self):
        x = [1.0, 0.0]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 1000 #Even with an enthalpy change, there should be no Gibbs free energy change for a pure substance
       
        expected_g = delta_h #Since there is no mixing, the Gibbs free energy change should be equal to the enthalpy change
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_delta_g_mix_with_delta_h_zero_enthalpy(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 0 #No enthalpy change, so it should reduce to the ideal case

        S = calculate_delta_s_mix(x, n=n, R=R)
        expected_g = -T * S
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_delta_g_mix_with_delta_h_negative_enthalpy(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = -1000 #Negative enthalpy change for mixing

        S = calculate_delta_s_mix(x, n=n, R=R)
        expected_g = delta_h - T * S
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_delta_g_mix_with_very_low_enthalpy(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n=1.0
        T = 300
        delta_h = -1e6 #Very low enthalpy change for mixing

        S = calculate_delta_s_mix(x, n=n, R=R)
        expected_g = delta_h - T * S
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_delta_g_mix_with_very_high_enthalpy(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 1e6 #Very high enthalpy change for mixing

        S = calculate_delta_s_mix(x, n=n, R=R)
        expected_g = delta_h - T * S
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_negative_temperature_with_delta_h(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = -50 #Negative temperature is not valid
        delta_h = 1000

        with self.assertRaises(ValueError):
            calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)

    def test_zero_temperature_with_delta_h(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 0 #At absolute zero, the Gibbs free energy change should be equal to the enthalpy change
        delta_h = 1000

        expected_g = delta_h
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_very_low_temperature_with_delta_h(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e-6
        delta_h = 1000

        expected_g = delta_h - T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_very_high_temperature_with_delta_h(self):
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e6
        delta_h = 1000

        expected_g = delta_h - T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

class TestOtherCases(unittest.TestCase):
    def test_invalid_gas_constant(self):
        x = [0.5, 0.5]
        R = 0 #Invalid gas constant
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)
    
    def test_invalid_type(self):
        x = "not a list" #Invalid type for mole fractions
        n = 1.0
        R = 8.314462618

        with self.assertRaises(TypeError):
            calculate_delta_s_mix(x, n=n, R=R)

if __name__ == '__main__':
    unittest.main()

#S = Entropy
#G = Gibbs free energy
#H = Enthalpy
#Delta S = R * sum(x_i * ln(x_i)) for i in components
#Delta G = Delta H - T * Delta S
#Delta H = 0 for ideal mixing, so Delta G = -T * Delta S
#Otherwise, Delta G = h - T * Delta S where h is the enthalpy change for mixing