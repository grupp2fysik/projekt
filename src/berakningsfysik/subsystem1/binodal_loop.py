"""En loop för att iterera över G för alla olika T,
och hitta kompositioner till binodalkurvan."""

import pandas as pd

df = pd.DataFrame(
    [[1, 2, 3], [4, 5, 3], [7, 8, 6], [9, 3, 10]],
    index=[0.2, 0.4, 0.6, 0.8],
    columns=["GT1", "GT2", "GT3"],
)

# använd df.iloc[:, start:slut] för att plocka ut G-kolumnerna ur
# den stora dataframen från termoberäkningarna

# för varje T: hitta konvext hölje - hitta fasett/tangent - hitta motsvarande x_a x_b -
# lägg x_a x_b i en array
for G in df.columns:
    data = df.loc[:, G]
    