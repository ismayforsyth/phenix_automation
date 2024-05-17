#!/dls/science/groups/i23/pyenvs/tihana_conda/bin/python

#set tab size to 2

import os
import re
from iotbx import reflection_file_reader 
import iotbx
import gemmi
import subprocess 
import copy 
from multiprocessing import Pool
from tqdm import tqdm
from halo import Halo # or can use tqdm, halo is prettier but don't tell you how long it will take
import plotly.graph_objects as go
import plotly.io as pio

cpus = os.cpu_count() - 1
pool = Pool(cpus)

os.environ['PHENIX'] = '/dls_sw/apps/phenix/1.20.1/phenix-1.20.1-4487'

mtzIn = input("File location for MTZ: ")
pdbIn = input("File location for PDB: ")
#seqIn = input("File location for SEQ: ")
projIn = input("Name of project: ")
# effIn = input("File location for EFF: ")
# pdbIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Phenix4/phaser_1/Lysozyme-FinalLSvsnLS_phaser.1.pdb"
# seqIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Lysozyme.seq"
# projIn = "Lysozyme"
# mtzIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/New_data_Clonly/LS/8keV/DataFiles/AUTOMATIC_DEFAULT_free.mtz"

elementsToTry = input("Which elements to try, comma separated: ")
#elementsToTry = "K, Cl, Ca, Xe"
elements = [x.strip() for x in elementsToTry.split(',')]

# Create mtz object storing information
mtzInfo = reflection_file_reader.any_reflection_file(mtzIn)
millerArray = mtzInfo.as_miller_arrays()
mtzInfo.file_content()
mtzobj = iotbx.mtz.object(mtzIn)

# Extract space group
space_group = mtzobj.space_group_name()
# Create dictionary matching concatenated space groups with correct spacing ???

# Extract unit cell params
csym = mtzobj.crystals()[0].crystal_symmetry()
unit_cell = csym.unit_cell()
# Remove list characters
unit_cell_strip = (str(unit_cell)).strip('()').replace(',', '')

# Extract WV
mtz = gemmi.read_mtz_file(mtzIn)
WV = (mtz.dataset(1).wavelength)
  
def wavelength_to_eV(wavelength):
    h = 6.62607015e-34
    c = 2.99792458e8
    eV = 1.602176634e-19
    energyeV = (h * c) / (float(wavelength) * 1e-10) / eV
    return energyeV

def lookup_fprime(energy, element):
    lookupPath = os.path.join("lookup", f"{element}.dat")
    closestEnergy = None
    closestValues = None
    with open(lookupPath, 'r') as lookupFile:
        for line in lookupFile:
            parts = line.split()
            currentEnergy, fPrime, fDoublePrime = map(float, parts)
            if closestEnergy is None or abs(energy - currentEnergy) < abs(energy - closestEnergy):
                closestEnergy = currentEnergy
                closestValues = (fPrime, fDoublePrime)
    return closestValues

def change_elem(pdbPath, elements):
	pdbInBase, pdbInExt = pdbPath.rsplit('.', 1)
	pdbInBase = os.path.basename(pdbInBase)
	print(f"Running on {str(pdbInBase)}")
	pdbLinesWrite = []
	pdbList = []

	with open(pdbPath, 'r') as file:
		pdbLinesWrite = file.readlines()

	toChange = []
	for line in pdbLinesWrite:
		if line.startswith("HETATM"):
			print(line.strip())
			changeHETATM = input("Would you like to run refinement on this HETATM? ")
			if changeHETATM.lower() in ('y', 'yes'):
				toChange.append(line)

	for element in elements:
		toFDPRefine = []
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
				toFDPRefine.append([element, line[21:22].strip(), line[23:26].strip()])
				
		with open(f"{pdbInBase}_{element}.{pdbInExt}", 'w') as pdbOut:
			for line in toWrite:
				if not line.endswith('\n'):
					line += '\n'
				pdbOut.write(line)
			pdbList.append((f"{pdbInBase}_{element}.{pdbInExt}", element, toFDPRefine))

	return pdbList

def run_phenix_refine(effFile):
	subprocess.run(["phenix.refine", effFile, "2>&1", "/dev/null"]) 
 
def scrapeLastAnomalousGroupData(log_file_path, theoreticalFDP):
  with open(log_file_path, 'r') as file:
    content = file.read()

  last_macro_cycle_index = content.rfind('MACRO_CYCLE')

  content_after_macro_cycle = content[last_macro_cycle_index:]
    
  pattern = re.compile(
      r'Anomalous scatterer group:\s+'
      r'Selection: "chain (?P<chain>[A-Za-z]) and resid (?P<resid>\d+) and element (?P<element>[A-Za-z]+)"\s+'
      r'Number of selected scatterers: \d+\s+'
      r'f_prime: +[+-]?\d+\.\d+\s+'
      r'f_double_prime: (?P<f_double_prime>[+-]?\d+\.\d+)'
  )

  matches = pattern.findall(content_after_macro_cycle)
  data = [(m[0], int(m[1]), m[2], float(m[3]), float(theoreticalFDP), float(abs(float(theoreticalFDP) - (float(m[3]))))) for m in matches]

  return data

def makeTable(scrapedData):
  header = ['Chain', 'ResidID', 'Element', 'Refined FDP', "Theoretical FDP", "Absolute Difference"]
  flattenedData = [item for sublist in scrapedData for item in sublist]
  tableData = list(map(list, zip(*flattenedData)))
  
  fig = go.Figure(data=[go.Table(
    header=dict(values=header,
                fill_color='paleturquoise',
                align='left'),
    cells=dict(values=tableData,
               fill_color='lavender',
               align='left'))
	])
  pio.write_html(fig, f'{projIn}.html')


# Create bpos eff file
def runBPos(pdbIn, elementIn):
  with open(f'bposEffParam_{elementIn}.eff', 'w') as file:
    file.write(f'''refinement {{
      crystal_symmetry {{  
        unit_cell = {unit_cell_strip}
        space_group = {space_group}
      }}
      input {{  
        pdb {{  
          file_name = "{pdbIn}"  
        }}  
        xray_data {{  
          file_name = "{mtzIn}"
          labels = IMEAN,SIGIMEAN  
          r_free_flags {{  
            file_name = "{mtzIn}"
            label = FreeR_flag  
            test_flag_value = 0  
          }}  
        }}  
        }}  
      }}  
      output {{  
        prefix = """{projIn}_bpos_{str(elementIn)}"""   
        job_title = """{projIn}"""  
        serial_format = "%d"
        write_def_file = False  
      }}  
      electron_density_maps {{  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT  
          mtz_label_phases = PH2FOFCWT  
          fill_missing_f_obs = True  
        }}  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT_no_fill  
          mtz_label_phases = PH2FOFCWT_no_fill  
        }}  
        map_coefficients {{  
          map_type = mFo-DFc  
          mtz_label_amplitudes = FOFCWT  
          mtz_label_phases = PHFOFCWT  
        }}  
        map_coefficients {{  
          map_type = anomalous  
          mtz_label_amplitudes = ANOM  
          mtz_label_phases = PHANOM  
        }}  
      }}  
      refine {{  
        strategy = *individual_sites individual_sites_real_space rigid_body \  
                  *individual_adp group_adp tls occupancies group_anomalous  
      }}  
      main {{  
        number_of_macro_cycles = 5  
        wavelength = {WV}
        nproc = {cpus}
      }}   
    }}''')

  subprocess.run(["phenix.refine", f"bposEffParam_{elementIn}.eff"])

def runFdp(elementIn, fPrime, toFDPRefine):
  anomalousScatterersStr = "anomalous_scatterers {\n"
  for elemID, chainID, residID in toFDPRefine:
      groupStr = f"""          group {{
            selection = chain {chainID} and resid {residID} and element {elemID}
            f_prime = {fPrime}
            refine = f_prime *f_double_prime
          }}\n"""
      anomalousScatterersStr += groupStr
  anomalousScatterersStr += "        }"
  
  with open(f'fdpEffParam_{elementIn}.eff', 'w') as file:
    file.write(f'''refinement {{  
			crystal_symmetry {{
        unit_cell = {unit_cell_strip}
        space_group = {space_group}
      }}  
      input {{  
        pdb {{  
          file_name = "{projIn}_bpos_{elementIn}_1.pdb"  
        }}  
        xray_data {{  
          file_name = "{mtzIn}"
          labels = I(+),SIGI(+),I(-),SIGI(-)
          r_free_flags {{  
            file_name = "{mtzIn}"
            label = FreeR_flag  
            test_flag_value = 0  
          }}  
        }}  
      }}  
      output {{  
        prefix = """{projIn}_fdp_{elementIn}"""  
        job_title = """{projIn}_{elementIn}"""  
        serial_format = "%d"
        write_def_file = False  
      }}  
      electron_density_maps {{  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT  
          mtz_label_phases = PH2FOFCWT  
          fill_missing_f_obs = True  
        }}  
        map_coefficients {{  
          map_type = 2mFo-DFc  
          mtz_label_amplitudes = 2FOFCWT_no_fill  
          mtz_label_phases = PH2FOFCWT_no_fill  
        }}  
        map_coefficients {{  
          map_type = mFo-DFc  
          mtz_label_amplitudes = FOFCWT  
          mtz_label_phases = PHFOFCWT  
        }}  
        map_coefficients {{  
          map_type = anomalous  
          mtz_label_amplitudes = ANOM  
          mtz_label_phases = PHANOM  
        }}  
      }}  
      refine {{  
        strategy = individual_sites individual_sites_real_space rigid_body \   
                  individual_adp group_adp tls occupancies *group_anomalous  
      	{anomalousScatterersStr}
      }}
      main {{  
        number_of_macro_cycles = 5  
        wavelength = {WV}
      }}  
    }}''')
	
  subprocess.run(["phenix.refine", f"fdpEffParam_{elementIn}.eff"])  

if __name__ == "__main__":    
	pdbList = change_elem(pdbIn, elements)
	energy = wavelength_to_eV(WV)
	scrapedData = []
	for pdb, ele, tfdpr in pdbList:
		print(pdb, ele, tfdpr)
		fp = list(lookup_fprime(energy, ele))[0]
		theoreticalFDP = list(lookup_fprime(energy, ele))[1]
		runBPos(pdb, ele)
		runFdp(ele, fp, tfdpr)
		logFile = (f'{projIn}_fdp_{ele}_1.log')
		scrapedData.append(scrapeLastAnomalousGroupData(logFile, theoreticalFDP))
	print(scrapedData)
	makeTable(scrapedData)

	






# with Halo(text="\nRunning B pos refinement", spinner="dots"):
#     run_phenix_refine("bposEffParam.eff")

# fDPRunList = []


# with Halo("\nRunning f'' refinement", spinner="dots"):
# 	pool.starmap(run_phenix_refine, fDPRunList)
