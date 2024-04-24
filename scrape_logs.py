import re

def scrapeLastAnomalousGroupData(log_file_path):
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
    data = [(m[0], int(m[1]), m[2], float(m[3])) for m in matches]

    return data

log_file_path = '/Users/vwg85559/phenix_automation-1/lys_fdp_refine_refine_2.log'
last_anomalous_group_data = scrapeLastAnomalousGroupData(log_file_path)
print(last_anomalous_group_data)

