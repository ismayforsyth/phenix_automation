import os
from iotbx import reflection_file_reader 
import cctbx
import iotbx
import gemmi
import subprocess 
import copy 
from multiprocessing import Pool
from tqdm import tqdm
from halo import Halo # or can use tqdm, halo is prettier but don't tell you how long it will take

cpus = os.cpu_count() - 1
pool = Pool(cpus)

os.system("module load phenix") # at some point need to do this nicer, eg. setenv()

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
			pdbList.append((pdbOut,element))
			for line in toWrite:
				if not line.endswith('\n'):
					line += '\n'
				pdbOut.write(line) 
    
	return pdbList

def new_func(pdbList):
    return pdbList 

def run_phenix_refine(effFile):
	subprocess.run(["phenix.refine", effFile, "2>&1", "/dev/null"]) 
   
# pdbIn = input("File location for PDB: ")
# seqIn = input("File location for SEQ: ")
# projIn = input("Name of project: ")
# effIn = input("File location for EFF: ")
pdbIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Phenix4/phaser_1/Lysozyme-FinalLSvsnLS_phaser.1.pdb"
seqIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Lysozyme.seq"
projIn = "Lysozyme"
# elementsToTry = input("Which elements to try, comma separated: ")
elementsToTry = "K, Cl, Ca, Xe"
elements = [x.strip() for x in elementsToTry.split(',')]


# mtzIn = input("File location for MTZ: ")
mtzIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/New_data_Clonly/LS/8keV/DataFiles/AUTOMATIC_DEFAULT_free.mtz"

# Create mtz object storing information
mtzInfo = reflection_file_reader.any_reflection_file(mtzIn)
millerArray = mtzInfo.as_miller_arrays()
mtzInfo.file_content()
mtzobj = iotbx.mtz.object(mtzIn)

# Extract space group
space_group = mtzobj.space_group_name()
print(space_group)
# Create dictionary matching concatenated space groups with correct spacing ???

# Extract unit cell params
csym = mtzobj.crystals()[0].crystal_symmetry()
unit_cell = csym.unit_cell()
# Remove list characters
unit_cell_strip = (str(unit_cell)).strip('()').replace(',', '')
print(unit_cell_strip)

# Extract WV
mtz = gemmi.read_mtz_file(mtzIn)
WV = (mtz.dataset(1).wavelength)
print(WV)

# Create bpos eff file
def runBPos(pdbIn,elementIn):
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
      sequence {{  
        file_name = "{seqIn}"  
      }}  
    }}  
    output {{  
      prefix = """{projIn}_bpos_{elementIn}"""   
      job_title = """{projIn}"""  
      serial_format = ""
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

def runFdp(elementIn):
  with open(f'fdpEffParam_{elementIn}.eff', 'w') as file:
    file.write(f'''refinement {{  
    crystal_symmetry {{
      unit_cell = {unit_cell_strip}
      space_group = {space_group}
    }}  
    input {{  
      pdb {{  
        file_name = "{projIn}_bpos_.pdb"  
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
      sequence {{  
        file_name = "{seqIn}"  
      }}  
    }}  
    output {{  
      prefix = """{projIn}_{elementIn}_fdp"""  
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
      anomalous_scatterers {{
        group {{
          selection = element {elementIn}
          f_prime = {fPrime}
          refine = f_prime *f_double_prime
        }}
      }}
    }}
    main {{  
      number_of_macro_cycles = 5  
      wavelength = {WV}
    }}  
  }}''')
	
  subprocess.run(["phenix.refine", f"fdpEffParam_{elementIn}.eff"])  

if __name__ == "__main__":
    pdbList = change_elem(pdbIn, elements)
    for pdb, ele in pdbList:
      runBPos(pdb, ele)
      runFdp(ele)

# with Halo(text="\nRunning B pos refinement", spinner="dots"):
#     run_phenix_refine("bposEffParam.eff")

# fDPRunList = []

# for element in elements:
#   energy = wavelength_to_eV(WV)
#   fPrime, fDoublePrime = lookup_fprime(energy, element)
#   fDPRunList.append(f"fdpEffParam_{element}.eff")

# with Halo("\nRunning f'' refinement", spinner="dots"):
# 	pool.starmap(run_phenix_refine, fDPRunList)
