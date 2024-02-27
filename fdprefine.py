import os
from iotbx import reflection_file_reader 
import cctbx
import iotbx
import gemmi
import subprocess 
import copy 

os.system("module load phenix")

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
			for line in toWrite:
				if not line.endswith('\n'):
					line += '\n'
				pdbOut.write(line)   
   
# pdbIn = input("File location for PDB: ")
# seqIn = input("File location for SEQ: ")
# projIn = input("Name of project: ")
# effIn = input("File location for EFF: ")
pdbIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Phenix4/phaser_1/Lysozyme-FinalLSvsnLS_phaser.1.pdb"
seqIn = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Lysozyme.seq"
projIn = "Lysozyme"

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
with open('/dls/i23/data/2023/cm33851-4/processing/Ismay/Scripts/bposEffParam.eff', 'w') as file:
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
    prefix = """{projIn}_bpos"""   
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
  }}   
}}''')

subprocess.run(["phenix.refine", "bposEffParam.eff"])
  
elements = ["Cl", "Ca", "K", "Mg"]

for element in elements:
  energy = wavelength_to_eV(WV)
  fPrime, fDoublePrime = lookup_fprime(energy, element)
  with open(f'fdpEffParam_{element}.eff', 'w') as file:
    file.write(f'''refinement {{  
    crystal_symmetry {{
      unit_cell = {unit_cell_strip}
      space_group = {space_group}
    }}  
    input {{  
      pdb {{  
        file_name = "{projIn}_bpos_1.pdb"  
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
      prefix = """{projIn}_fdp"""  
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
      strategy = individual_sites individual_sites_real_space rigid_body \   
                individual_adp group_adp tls occupancies *group_anomalous  
      anomalous_scatterers {{
        group {{
          selection = element {element}
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

subprocess.run(["phenix.refine", "fdpEffParam.eff"])  

# declare residues and give theoretical f' etc --- crossec (f'/f'' tabulated values)

# run f'' using phenix.refine

# subprocess.run(["phenix.refine", "fdpEffParam.eff"])

# change element type in pdb and redo.... 









# os.system("module load phenix")
# os.popen("phenix", "bposEffParam.eff")

# cwd = os.getcwd()


# # phenix refine input setup
# effAnomParams = input("What Phenix group selections do you want to change and refine? Comma separated. ")
# effAnomParamsList = [x.strip() for x in effAnomParams.split(',')]

# # pdb file input setup
# pdbSites = input("What are the sites called in the PDB that you want to change and refine? Comma separated. ")
# pdbSitesList = [x.strip() for x in pdbSites.split(',')]

# elementsToTry = input("Which elements to try, comma separated: ")
# elementsToTryList = [x.strip() for x in elementsToTry.split(',')]

# # modify eff file
# for element in elementsToTryList:
#     effOutputFile = f"refine_{element}.eff"
#     with open(effIn, 'r') as effInFile, open(os.path.join(cwd+ effOutputFile), 'w') as outfile:
#         for line in effInFile:
#             if line.startswith("        selection = chain"):
#                 for term in effAnomParamsList:
#                     if line.endswith(term):
#                         pass
#                         # now change line below to whatever f
#                     else:
#                         pass
#             elif line.startswith("        selection = element"):
#                 for term in effAnomParamsList:
#                     if line.endswith(term):
#                         pass
#                         # changle line below
#                         # write new element in this line
#                     else:
#                         pass
#             else:
#                 outfile.write(line)



# refinement {
#   crystal_symmetry {
#     unit_cell = 79.1429 79.1429 36.9147 90 90 90
#     space_group = P 43 21 2
#   }
#   input {
#     pdb {
#       file_name = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Phenix4/phaser_1/Lysozyme-FinalLSvsnLS_phaser.1.pdb"
#     }
#     xray_data {
#       file_name = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/New_data_Clonly/LS/8keV/DataFiles/AUTOMATIC_DEFAULT_free.mtz"
#       labels = IMEAN,SIGIMEAN
#       r_free_flags {
#         file_name = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/New_data_Clonly/LS/8keV/DataFiles/AUTOMATIC_DEFAULT_free.mtz"
#         label = FreeR_flag
#         test_flag_value = 0
#       }
#     }
#     sequence {
#       file_name = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Lysozyme.seq"
#     }
#   }
#   output {
#     prefix = """Lysozyme-FinalLSvsnLS_refine"""
#     serial = 4
#     serial_format = "%d"
#     job_title = """Lysozyme LS BPos - 8keV"""
#     write_def_file = False
#   }
#   electron_density_maps {
#     map_coefficients {
#       map_type = 2mFo-DFc
#       mtz_label_amplitudes = 2FOFCWT
#       mtz_label_phases = PH2FOFCWT
#       fill_missing_f_obs = True
#     }
#     map_coefficients {
#       map_type = 2mFo-DFc
#       mtz_label_amplitudes = 2FOFCWT_no_fill
#       mtz_label_phases = PH2FOFCWT_no_fill
#     }
#     map_coefficients {
#       map_type = mFo-DFc
#       mtz_label_amplitudes = FOFCWT
#       mtz_label_phases = PHFOFCWT
#     }
#     map_coefficients {
#       map_type = anomalous
#       mtz_label_amplitudes = ANOM
#       mtz_label_phases = PHANOM
#     }
#   }
#   refine {
#     strategy = *individual_sites individual_sites_real_space rigid_body 
#                *individual_adp group_adp tls occupancies group_anomalous
#   }
#   main {
#     number_of_macro_cycles = 5
#     wavelength = 1.5498
#   }
#   gui {
#     base_output_dir = /dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Phenix4
#     tmp_dir = "/dls/i23/data/2023/cm33851-4/processing/Ismay/Lysozyme/Phenix4/.phenix/tmp"
#     notify_email = """ismay.forsyth@diamond.ac.uk"""
#   }
# }
