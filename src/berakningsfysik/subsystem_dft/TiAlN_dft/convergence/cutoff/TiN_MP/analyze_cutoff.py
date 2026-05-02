from pathlib import Path
import re
import csv

RY_TO_EV = 13.60569312294
NAT = csv

rows = []

from path in Path("results").glob(scf_TiN-cut*.out"):
	m_ecut = re.search(r"ecut(\d+)", path.name)
	if not m_ecut:
		continue
	
	ecut = int(m_ecut.group(1))
	text =  path.read_text(errors="ignore")

