import os
from iotbx import reflection_file_reader 
import cctbx
import iotbx
import gemmi
import subprocess 

os.system("module load phenix")

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
  
# create fdp eff file
with open('fdpEffParam.eff', 'w') as file:
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
