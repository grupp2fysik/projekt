"""Denna fil innehåller enhetstester för funktionerna calculate_delta_s_mix, calculate_delta_g_mix
och calculate_delta_g_mix_with_h i thermo.py. Testerna täcker olika scenarier, 
inklusive ideal blandning, olika temperaturer, olika antal mol, ogiltiga indata och mer. 

Syftet är att säkerställa att funktionerna fungerar korrekt under olika 
förhållanden och hanterar ogiltiga indata på ett lämpligt sätt."""

import unittest
import math
from thermo import calculate_delta_s_mix, calculate_delta_g_mix, calculate_delta_g_mix_with_h

class TestDeltaSMix(unittest.TestCase):
    """
    En test för klassen TestDeltaSMix som innehåller enhetstester för funktionen calculate_delta_s_mix.
    Testerna inkluderar olika scenarier för att säkerställa att funktionen fungerar korrekt
    under olika förhållanden och hanterar ogiltiga indata på ett lämpligt sätt.
    """
    def test_ideal_mixing(self):
        """
        Testa att calculate_delta_s_mix returnerar det förväntade värdet för en 
        ideal blandning av två komponenter i lika proportioner.

        Parametrar:
        x (list): En lista med molbråk för komponenterna i blandningen.
        R (float): Gas konstanten (standard är 8.314462618 J/(mol*K)).
        n (float): Det totala antalet mol i blandningen.
        """

        x = [0.5, 0.5] 
        R = 8.314462618 
        n = 1.0 

        expected_s = -R * (0.5 * math.log(0.5) + 0.5 * math.log(0.5)) #Förväntad entropi av blandning för en ideal blandning av två komponenter i lika proportioner.
        calculated_s = calculate_delta_s_mix(x, n=n, R=R) #Säkerställa att den beräknade entropin av blandning är lika med det förväntade värdet, med en hög precision på 100 decimaler.
        self.assertAlmostEqual(calculated_s, expected_s, places=100) #places=100 säkerställer att jämförelsen görs med en precision på 100 decimaler, vilket är viktigt för att verifiera noggrannheten i beräkningen av entropi av blandning.


    def test_delta_s_bigger_than_zero(self):
        """
        Testa att calculate_delta_s_mix returnerar ett värde större än noll för en blandning av två komponenter i lika proportioner,
        vilket är förväntat eftersom blandning av två olika komponenter bör öka oordningen och därmed entropin.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0

        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertGreater(calculated_s, 0) #Assert that the calculated entropy change is greater than zero, which is expected for a mixing process that increases disorder.

    def test_different_mole_count(self):
        """
        Testa ett fall där det totala antalet mol i blandningen är större än 1,
        vilket bör resultera i en proportionellt större entropi av blandning.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 2.0

        expected_s = -R * n * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

    def test_more_than_one_component(self):
        """
        Testa ett fall där det finns mer än två komponenter i blandningen,
        vilket bör resultera i en entropi av blandning som tar hänsyn till alla komponenter.
        """

        x = [0.25, 0.25, 0.25, 0.25]
        R = 8.314462618
        n = 1.0

        expected_s = -R * 4 * (0.25 * math.log(0.25)) #Each component contributes equally to the entropy
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)
    
    def test_division_of_three(self):
        """
        Testa ett fall där det finns tre komponenter i blandningen i lika proportioner,
        vilket bör inte gå för att summan av molbråken inte är lika med 1, och därmed bör resultera i ett ValueError.
        """

        x = [0.333333333333333, 0.333333333333333, 0.333333333333333]
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #Sum != 1

    def test_pure_substance(self):
        """
        Testa ett fall där det finns en ren substans i blandningen, vilket bör resultera i en entropi av blandning lika med noll,
        eftersom det inte finns någon blandning och därmed ingen ökning av oordningen.
        """

        x = [1.0, 0.0]
        R = 8.314462618
        n = 1.0

        self.assertEqual(calculate_delta_s_mix(x, n=n, R=R), 0.0)
        #Equal to 0 because there is no mixing, only one component
    
    def test_mole_fractions_not_equal_to_one(self):
        """
        Testa ett fall där summan av molbråken inte är lika med 1,
        vilket bör resultera i ett ValueError.
        """

        x = [0.5, 0.2] #Sum != 1
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #Sum != 1

    def test_negative_mole_fractions(self):
        """
        Testa ett fall där det finns negativa molbråken i blandningen,
        vilket bör resultera i ett ValueError.
        """

        x = [-0.1, 1.1] #There is a negative mole fraction, which is not valid
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)
    
    def test_invalid_mole_fractions(self):
        """
        Testa ett fall där det finns mer än två komponenter i blandningen i lika proportioner,
        vilket bör inte gå för att summan av molbråken inte är lika med 1, och därmed bör resultera i ett ValueError.
        """

        x = [0.5, 0.5, 0.5] #Sum > 1
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R) #Sum != 1

    def test_no_mole_fractions(self):
        """
        Testa ett fall där det inte finns några molbråken i blandningen,
        vilket bör resultera i ett ValueError.
        """

        x = []
        R = 8.314462618
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)
    
    def test_almost_one(self):
        """
        Testa ett fall där det finns en komponent i blandningen som är nästan lika med 1, 
        och den andra komponenten är nästan lika med 0, vilket bör resultera i en entropi av blandning som är mycket
        """

        x = [0.9999999999999999, 1e-16] #Sum = 1, men en av komponenterna är obetydlig.
        R = 8.314462618
        n = 1.0

        expected_s = -R * (0.9999999999999999 * math.log(0.9999999999999999) + 1e-16 * math.log(1e-16))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)
 
    def test_negative_n(self):
        """
        Testa ett fall där molen i blandningen är negativt, 
        vilket bör resultera i ett ValueError eftersom det inte är fysiskt meningsfullt att ha ett negativt antal mol i en blandning.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = -1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)

    def test_zero_n(self):
        """
        Testa ett fall där molen i blandningen är noll,
        vilket bör resultera i en entropi av blandning som är noll,
        eftersom det inte finns någon blandning och därmed ingen ökning av oordningen.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 0

        expected_s = 0.0
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

    def test_very_high_n(self):
        """
        Testa ett fall där molen i blandningen är mycket hög,
        vilket bör resultera i en proportionellt större entropi av blandning,
        eftersom entropin av blandning är proportionell mot antalet mol i blandningen.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1e6

        expected_s = -R * n * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

    def test_very_low_n(self):
        """
        Testa ett fall där molen i blandningen är mycket låg,
        vilket bör resultera i en proportionellt mindre entropi av blandning,
        eftersom entropin av blandning är proportionell mot antalet mol i blandningen.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1e-6 #Very low n means a very small entropy change

        expected_s = -R * n * (0.5 * math.log(0.5) + 0.5 * math.log(0.5))
        calculated_s = calculate_delta_s_mix(x, n=n, R=R)
        self.assertAlmostEqual(calculated_s, expected_s, places=100)

class TestDeltaGMix(unittest.TestCase):
    """
    En test för klassen TestDeltaGMix som innehåller enhetstester för funktionen calculate_delta_g_mix.
    Testerna inkluderar olika scenarier för att säkerställa att funktionen fungerar korrekt
    under olika förhållanden och hanterar ogiltiga indata på ett lämpligt sätt.
    """

    def test_delta_g_mix_temperature(self):
        """
        Testa att calculate_delta_g_mix returnerar det förväntade värdet
        för en ideal blandning av två komponenter i lika proportioner vid en given temperatur.

        Parametrar:
        x (list): En lista med molbråk för komponenterna i blandningen.
        R (float): Gas konstanten (standard är 8.314462618 J/(mol*K)).
        n (float): Det totala antalet mol i blandningen.
        T (float): Temperaturen i Kelvin.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300 #Tempature in Kelvin

        expected_g = -T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_pure_substance(self):
        """
        Testa ett fall där det finns en ren substans i blandningen, vilket bör resultera i en Gibbs fri energi av blandning lika med noll,
        eftersom det inte finns någon blandning och därmed ingen förändring i Gibbs fri energi.
        """

        x = [1.0, 0.0]
        R = 8.314462618
        n = 1.0
        T = 300

        expected_g = 0.0 #Ingen förändring i Gibbs fri energi eftersom det inte finns någon blandning, bara en komponent
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_negative_temperature(self):
        """
        Testa ett fall där temperaturen är negativ, vilket bör resultera i ett ValueError 
        eftersom det inte är fysiskt meningsfullt att ha en negativ temperatur i Kelvin.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        T = -50
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_g_mix(x, T, n=n, R=R)

    def test_zero_temperature(self):
        """
        Testa ett fall där temperaturen är noll, vilket bör resultera i en Gibbs fri energi av blandning lika med noll,
        eftersom det inte finns någon förändring i Gibbs fri energi vid absolut nollpunkt.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 0

        expected_g = 0.0
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_very_low_temperature(self):
        """
        Testa ett fall där temperaturen är mycket låg, vilket bör resultera i en mycket liten Gibbs fri energi av blandning,
        eftersom Gibbs fri energi av blandning är proportionell mot temperaturen.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e-6

        expected_g = -T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_very_high_temperature(self):
        """
        Testa ett fall där temperaturen är mycket hög, vilket bör resultera i en mycket stor Gibbs fri energi av blandning,
        eftersom Gibbs fri energi av blandning är proportionell mot temperaturen.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e6

        expected_g = -T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix(x, T, n=n, R=R)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_temperature_relationship(self):
        """
        Testa att calculate_delta_g_mix returnerar ett värde som är proportionellt mot temperaturen,
        vilket är förväntat eftersom Gibbs fri energi av blandning är proportionell mot temperaturen.

        Parametrar:
        T1 (float): En temperatur i Kelvin.
        T2 (float): En annan temperatur i Kelvin som är dubbelt så hög som T
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T1 = 300
        T2 = 600

        g1 = calculate_delta_g_mix(x, T1, n=n, R=R)
        g2 = calculate_delta_g_mix(x, T2, n=n, R=R)
        self.assertAlmostEqual(g2, 2 * g1, places=100)

class TestDeltaGMixWithH(unittest.TestCase):
    """
    En test för klassen TestDeltaGMixWithH som innehåller enhetstester för funktionen calculate_delta_g_mix_with_h.
    """

    def test_delta_g_mix_with_delta_h(self):
        """
        Testa att calculate_delta_g_mix_with_h returnerar det förväntade värdet
        för en ideal blandning av två komponenter i lika proportioner vid en given temperatur,
        när entalpi av blandning är känd.

        Parametrar:
        x (list): En lista med molbråk för komponenterna i blandningen.
        R (float): Gas konstanten (standard är 8.314462618 J/(mol*K)).
        n (float): Det totala antalet mol i blandningen.
        T (float): Temperaturen i Kelvin.
        delta_h (float): Entalpi av blandning i J/mol.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 1000

        expected_g = delta_h - T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_delta_g_mix_with_delta_h_pure_substance(self):
        """
        Testa ett fall där det finns en ren substans i blandningen, vilket bör resultera i en Gibbs fri energi av blandning lika med entalpi av blandning,
        eftersom det inte finns någon blandning och därmed ingen förändring i Gibbs fri energi, 
        men entalpi av blandning kan fortfarande vara närvarande på grund av interaktioner mellan molekylerna i den rena substansen.
        """

        x = [1.0, 0.0]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 1000
       
        expected_g = delta_h
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_delta_g_mix_with_delta_h_zero_enthalpy(self):
        """
        Testa att calculate_delta_g_mix_with_h returnerar det förväntade värdet
        för en ideal blandning av två komponenter i lika proportioner vid en given temperatur,
        när entalpi av blandning är noll.
        """
        
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = 0

        S = calculate_delta_s_mix(x, n=n, R=R)
        expected_g = -T * S
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_delta_g_mix_with_delta_h_negative_enthalpy(self):
        """
        Testa att calculate_delta_g_mix_with_h returnerar det förväntade värdet
        för en ideal blandning av två komponenter i lika proportioner vid en given temperatur,
        när entalpi av blandning är negativ.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 300
        delta_h = -1000

        S = calculate_delta_s_mix(x, n=n, R=R)
        expected_g = delta_h - T * S
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100) #Exotermisk blandning
    
    def test_delta_g_mix_with_very_low_enthalpy(self):
        """
        Testa att calculate_delta_g_mix_with_h returnerar det förväntade värdet
        för en ideal blandning av två komponenter i lika proportioner vid en given temperatur,
        när entalpi av blandning är mycket låg.
        """

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
        """
        Testa att calculate_delta_g_mix_with_h returnerar det förväntade värdet
        för en ideal blandning av två komponenter i lika proportioner vid en given temperatur,
        när entalpi av blandning är mycket hög.
        """

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
        """
        Testa ett fall där temperaturen är negativ, vilket bör resultera i ett ValueError
        eftersom det inte är fysiskt meningsfullt att ha en negativ temperatur i Kelvin.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = -50 #Negative temperature is not valid
        delta_h = 1000

        with self.assertRaises(ValueError):
            calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)

    def test_zero_temperature_with_delta_h(self):
        """
        Testa ett fall där temperaturen är noll, vilket bör resultera i en 
        Gibbs fri energi av blandning lika med entalpi av blandning.
        """
        
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 0
        delta_h = 1000

        expected_g = delta_h
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)
    
    def test_very_low_temperature_with_delta_h(self):
        """
        Testa ett fall där temperaturen är mycket låg, vilket bör resultera i en Gibbs fri energi av 
        blandning som är nära entalpi av blandning, eftersom den term som involverar entropi av blandning
        blir mycket liten vid låga temperaturer.
        """
        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e-6
        delta_h = 1000

        expected_g = delta_h - T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

    def test_very_high_temperature_with_delta_h(self):
        """
        Testa ett fall där temperaturen är mycket hög, vilket bör resultera i en Gibbs fri energi av 
        blandning som är nära entalpi av blandning, eftersom den term som involverar entropi av blandning
        blir mycket liten vid höga temperaturer.
        """

        x = [0.5, 0.5]
        R = 8.314462618
        n = 1.0
        T = 1e6
        delta_h = 1000

        expected_g = delta_h - T * calculate_delta_s_mix(x, n=n, R=R)
        calculated_g = calculate_delta_g_mix_with_h(x, T, n=n, R=R, delta_h=delta_h)
        self.assertAlmostEqual(calculated_g, expected_g, places=100)

class TestOtherCases(unittest.TestCase):
    """
    En test för klassen TestOtherCases som innehåller enhetstester för...
    andra saker som inte finns redan i de andra testklasserna.
    """

    def test_invalid_gas_constant(self):
        """
        Testa ett fall där gas konstanten är inte 8.314462618,
        vilket bör resultera i ett ValueError eftersom det inte är en giltig gas konstant.
        """

        x = [0.5, 0.5]
        R = 0
        n = 1.0

        with self.assertRaises(ValueError):
            calculate_delta_s_mix(x, n=n, R=R)
    
    def test_invalid_type(self):
        """
        Testa ett fall där molbråken inte är en lista, vilket bör resultera i ett TypeError
        eftersom molbråken måste vara en lista av flyttal som representerar molbråken för varje komponent i blandningen.
        """

        x = "not a list"
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
#Delta H = 0 för ideal mixing, så Delta G = -T * Delta S
#Annars, Delta G = h - T * Delta S där h är entalpi av blandning.
