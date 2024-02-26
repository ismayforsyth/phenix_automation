import copy

elements = ["Na", "Mg", "S", "Cl", "K", "Ca"]

pdbIn = "/Users/vwg85559/phenix_automation/Lysozyme-FinalLSvsnLS_phaser.1.pdb"
pdbInBase, pdbInExt = pdbIn.rsplit('.', 1)

pdbLinesWrite = []

with open(pdbIn, 'r') as file:
    pdbLinesWrite = file.readlines()
    pdbLinesWrite_const = pdbLinesWrite

toChange = []
for line in pdbLinesWrite:
    if line.startswith("HETATM"):
        print(line.strip())
        changeHETATM = input("Would you like to run refinement on this HETATM? ")
        if changeHETATM.lower() in ('y', 'yes'):
            toChange.append(line)
        
for element in elements:
    toWrite = copy.deepcopy(pdbLinesWrite)
    for i, line in enumerate(pdbLinesWrite):
        if line in toChange:
            elementSymbol = element.rjust(2)[:2] if len(element) == 2 else ' ' + element
            residueName = element.rjust(3)[:3]
            elementName = element.ljust(4)[:4]
            lineList = list(line)
            lineList[12:16] = elementName
            lineList[17:20] = residueName
            lineList[76:78] = elementSymbol
            toWrite[i] = ''.join(lineList)
            
    with open(f"{pdbInBase}_{element}.{pdbInExt}", 'w') as pdbOut:
        for line in toWrite:
            if not line.endswith('\n'):
                line += '\n'
            pdbOut.write(line)