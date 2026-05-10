--INTRODUKTION--
Detta är ett program för att skapa ett fasdiagram för en binär legering.
Programmet kan dessutom plotta blandingsentalpi, den konfigurationella blandningsentropin
samt gibbs fria blandningsenergi. Därutöver kan programmet räkna ut kompositionen och 
andelen av varje fas vid en given temperatur. Programmet använder utdata-filer från
DFT-beräkningar som gjort med mjukvaran Quantum Espresso.


--INNAN DU GÖR BERÄKNINAGR--
Se först till att QE-output-filer för din legering finns någonstans i mappen "Enthalpy_interpolation".
Innan några beräkningar kan göras måste programmet förses med en csv-fil som innehåller 
materialspecifika parametrar. Filen ska döpas till legeringens namn, till exempel "TiAlN.csv", 
och placeras i undermappen "alloy_parameters".

Nedan följer ett exempel på hur denna fil kan se ut:

atomer per metallplats, 2
temperatur min, 0
temperatur max, 10000
filväg till Quantum Espresso-filer, TiAlN
specifika temperaturer, 1000 5000

Här kommer en förklaring av parametrarna:

n = antal atomer per metallplats i legeringens kristallstruktur
temperatur min = lägsta temperatur som plottas i fasdiagrammet
temperatur max = högsta temperatur som plottas i fasdiagrammet
filväg till Quantum Espresso-filer = filväg, utgående från undermappen "enthalpy_interpolation", där output
                                     från Quantum Espresso finns
speciella temperaturer = om man vill få ut information, till exempel en plott av Gibbs fria 
                        blandningsenergi, om legeringen vid en eller flera 
                        specifika temperaturer kan man ange dessa här. 


--HUR DU GÖR BERÄKNINGAR--

1. Innan några annat kan köras, så ska följande kommando skrivas i teminalen:

>> ./calculate.sh <legering>

där <legering> är namnet på din legering, till exempel TiAlN.

Detta kommando kommer interpolera blandingsentalpin, beräkna de andra termodynamiska storheterna,
och sedan med hjälp av dessa hitta binodal- och spinodalkurvor för att plotta fasdiagrammet. Sedan plottas 
blandingsentalpin och blandningsentropin. Figurerna läggs i undermappen "plots".

2. För att ta reda på fasandelar vid en eller flera temperaturer körs följande kommando:

>> ./phase_analysis.sh <legering> <temp1> <temp2>

där <legering> är namnet på den legering, och <temp1> och <temp2> är de temperaturer du vill 
ska analyseras. Du kan inkludera hur många temperaturer du vill: antingen bara en eller flera.
Men det är viktigt att dessa temperaturer finns inlagda som "specifika temperaturer" i din parameter-fil
i mappen "alloy_paramters". 

Nu kommer en fasanalys att skrivas ut i terminalen för dina valda temperaturer. 