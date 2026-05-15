--INTRODUKTION--
Detta är ett program för att skapa ett fasdiagram för en binär legering.
Programmet kan dessutom plotta blandingsentalpi, den konfigurationella blandningsentropin
samt gibbs fria blandningsenergi. Därutöver kan programmet räkna ut kompositionen och 
andelen av varje fas vid en given temperatur. Programmet använder utdata-filer från
DFT-beräkningar som gjorts med mjukvaran Quantum Espresso (QE).

OBS: Om du inte laddat ner uv, så kör de versionerna av bash-scripten som slutar med "_py.sh".

--INNAN DU GÖR BERÄKNINAGR--
Se först till att en csv-fil med data för beräknad blandningsentalpi finns tillgänglig. Filen bör ha
en kolumn "x" för undersökta kompositioner, och en kolumn "hmix_ev_per_atom" för beräknad blandningsentalpi 
vid dessa kompositioner. 

Innan några beräkningar kan göras måste programmet förses med en csv-fil som innehåller 
materialspecifika parametrar. Filen ska döpas till legeringens namn, till exempel "TiAlN.csv", 
och placeras i undermappen "subsystem_python/alloy_parameters".

Nedan följer ett exempel på hur denna fil ska se ut:

atomer per metallplats, 2
temperatur min, 0
temperatur max, 10000
filväg till Quantum Espresso-filer, ../subsystem_dft/results/final_hmix_sqs/final_fixed_sqs_hmix_results.csv
specifika temperaturer, 1000 5000

Här kommer en förklaring av parametrarna:

n: Antal atomer per metallplats i legeringens kristallstruktur.

temperatur min: Lägsta temperatur som kommer analyseras (anges som heltal grader K).

temperatur max: Högsta temperatur som kommer analyseras (anges som heltal grader K).

filväg till h_mix-fil: Filväg till csv-fil med beräknad blandningsentalpi.

speciella temperaturer: Om man vill få ut information, till exempel en plott av Gibbs fria 
                        blandningsenergi, om legeringen vid en eller flera 
                        specifika temperaturer kan man ange dessa här. Dessa temperaturer bör ligga i 
                        intervallet [temperatur min, temperatur max] och anges som grader K.


--HUR DU GÖR BERÄKNINGAR--

1. Innan några annat kan köras, så ska följande kommando skrivas i teminalen:

>> ./calculate.sh <legering>

där <legering> är namnet på din legering, till exempel TiAlN.

Detta kommando kommer interpolera blandingsentalpin, beräkna de andra termodynamiska storheterna,
och sedan med hjälp av dessa hitta binodal- och spinodalkurvor för att plotta fasdiagrammet. Sedan plottas fasdiagram,
blandingsentalpin och blandningsentropin. Figurerna läggs i undermappen "plots".

2. För att ta reda på fasandelar vid en eller flera temperaturer körs följande kommando: (obs just nu görs analys för alla temperaturer)

>> ./phase_analysis.sh <legering> <temp1> <temp2>

där <legering> är namnet på den legering, och <temp1> och <temp2> är de temperaturer du vill 
ska analyseras. Du kan inkludera hur många temperaturer du vill: antingen bara en eller flera.
Men det är viktigt att dessa temperaturer fanns inlagda som "specifika temperaturer" i din parameter-fil
i mappen "alloy_paramters" under körningen av ./calculate.sh.

Nu kommer en fasanalys att skrivas ut i terminalen för dina valda temperaturer. 

3. För att plotta Gibbs fria blandningsenergi för en specifik temperatur körs följande kommando i terminalen:

>> ./plot_gibbs.sh <legering> <temp1> <temp2>

där <legering> är namnet på den legering, och <temp1> och <temp2> är de temperaturer du vill 
ska analyseras. Du kan inkludera hur många temperaturer du vill: antingen bara en eller flera.
Men det är viktigt att dessa temperaturer fanns inlagda som "specifika temperaturer" i din parameter-fil
i mappen "alloy_parameters" under körningen av ./calculate.sh.

Nu kommer en plot skapas för varje temperatur och läggas i undermappen "plots".