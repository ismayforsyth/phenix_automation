import plotly.graph_objects as go
import plotly.io as pio

data = [[('C', 1, 'Mg', -4.60658, 1.635705, 6.242285), ('C', 2, 'Mg', -0.601712, 1.635705, 2.2374169999999998), ('C', 3, 'Mg', -1.84085, 1.635705, 3.4765550000000003), ('C', 4, 'Mg', 5.39997, 1.635705, 3.764265), ('J', 1, 'Mg', -26.4326, 1.635705, 28.068305000000002), ('J', 2, 'Mg', 5.43078, 1.635705, 3.7950750000000006), ('J', 3, 'Mg', -20.4037, 1.635705, 22.039405000000002), ('J', 4, 'Mg', -39.1703, 1.635705, 40.806005), ('J', 5, 'Mg', -11.4743, 1.635705, 13.110005), ('J', 7, 'Mg', 1.29384, 1.635705, 0.34186499999999986), ('J', 9, 'Mg', -18.3873, 1.635705, 20.023005), ('D', 1, 'Mg', 0.722764, 1.635705, 0.912941), ('D', 2, 'Mg', 19.4999, 1.635705, 17.864195)], [('C', 1, 'P', 3.16924, 3.583364, 0.41412400000000016), ('C', 2, 'P', 2.39578, 3.583364, 1.1875840000000002), ('C', 3, 'P', -3.29834, 3.583364, 6.881704), ('C', 4, 'P', 7.1598, 3.583364, 3.5764359999999997), ('J', 1, 'P', -18.7704, 3.583364, 22.353763999999998), ('J', 2, 'P', 4.83474, 3.583364, 1.251376), ('J', 3, 'P', -15.8458, 3.583364, 19.429164), ('J', 4, 'P', -34.9369, 3.583364, 38.520264000000005), ('J', 5, 'P', -6.28578, 3.583364, 9.869144), ('J', 7, 'P', 0.207765, 3.583364, 3.375599), ('J', 9, 'P', -18.4161, 3.583364, 21.999464), ('D', 1, 'P', 5.2535, 3.583364, 1.6701359999999998), ('D', 2, 'P', 22.8577, 3.583364, 19.274336)], [('C', 1, 'S', 9.90908, 0.4049515, 9.5041285), ('C', 2, 'S', 8.86613, 0.4049515, 8.4611785), ('C', 3, 'S', -5.71518, 0.4049515, 6.1201315), ('C', 4, 'S', -2.32787, 0.4049515, 2.7328215), ('J', 1, 'S', -28.2012, 0.4049515, 28.6061515), ('J', 2, 'S', 6.76159, 0.4049515, 6.3566385), ('J', 3, 'S', -22.4151, 0.4049515, 22.820051499999998), ('J', 4, 'S', -45.0168, 0.4049515, 45.421751500000006), ('J', 5, 'S', -14.577, 0.4049515, 14.9819515), ('J', 7, 'S', 2.48114, 0.4049515, 2.0761884999999998), ('J', 9, 'S', -18.8407, 0.4049515, 19.245651499999997), ('D', 1, 'S', 5.23668, 0.4049515, 4.8317285), ('D', 2, 'S', 23.3354, 0.4049515, 22.9304485)], [('C', 1, 'K', 3.7068, 0.8579073, 2.8488927), ('C', 2, 'K', 5.67925, 0.8579073, 4.8213427), ('C', 3, 'K', -1.06324, 0.8579073, 1.9211473), ('C', 4, 'K', 9.46155, 0.8579073, 8.6036427), ('J', 1, 'K', -35.375, 0.8579073, 36.2329073), ('J', 2, 'K', 4.42682, 0.8579073, 3.5689127000000003), ('J', 3, 'K', -25.6455, 0.8579073, 26.5034073), ('J', 4, 'K', -57.1551, 0.8579073, 58.0130073), ('J', 5, 'K', -22.7475, 0.8579073, 23.6054073), ('J', 7, 'K', -3.22055, 0.8579073, 4.0784573), ('J', 9, 'K', -13.5062, 0.8579073, 14.3641073), ('D', 1, 'K', 7.24921, 0.8579073, 6.3913027), ('D', 2, 'K', 25.8151, 0.8579073, 24.9571927)], [('C', 1, 'Cl', 4.83882, 0.5274955, 4.3113245000000004), ('C', 2, 'Cl', -0.8546, 0.5274955, 1.3820955000000001), ('C', 3, 'Cl', -4.97009, 0.5274955, 5.4975855), ('C', 4, 'Cl', -1.56921, 0.5274955, 2.0967055), ('J', 1, 'Cl', -26.6076, 0.5274955, 27.135095500000002), ('J', 2, 'Cl', 2.46038, 0.5274955, 1.9328844999999997), ('J', 3, 'Cl', -24.607, 0.5274955, 25.1344955), ('J', 4, 'Cl', -43.8992, 0.5274955, 44.4266955), ('J', 5, 'Cl', -10.6853, 0.5274955, 11.2127955), ('J', 7, 'Cl', 5.66495, 0.5274955, 5.1374545000000005), ('J', 9, 'Cl', -20.9667, 0.5274955, 21.4941955), ('D', 1, 'Cl', 11.6594, 0.5274955, 11.1319045), ('D', 2, 'Cl', 24.9459, 0.5274955, 24.4184045)]]

import matplotlib.pyplot as plt
import pandas as pd

df = pd.DataFrame(data, columns=['C1', 'C2', 'C3', 'C4', 'J1', 'J2', 'J3', 'J4', 'J5', 'J7', 'J9', 'D1', 'D2'] )
print(df)

#['Chain', 'ResidID', 'Element', 'Refined FDP', "Theoretical FDP", "Absolute Difference"]
#'Mg', 'P', 'S', 'K', 'Cl'

all_data = []
for chain_id in data:
    id = pd.DataFrame(chain_id, columns=['Chain', 'ResidID', 'Element', 'Refined FDP', "Theoretical FDP", "Absolute Difference"])
    all_data.append(id)

# Concatenating all the dataframes
combined_data = pd.concat(all_data)
#combined_data['new_column'] = combined_data['Chain'] + str(combined_data['ResidID'])
print(combined_data['ResidID'])

fig = go.Figure(data=[go.Table(
    header=dict(values=list(combined_data.columns),
                fill_color='paleturquoise',
                align='left'),
    cells=dict(values=combined_data.transpose().values.tolist(),
               fill_color='lavender',
               align='left'))
])

fig.show()

# def makeTable(scrapedData):
#   header = ['Chain', 'ResidID', 'Element', 'Refined FDP', "Theoretical FDP", "Absolute Difference"]
#   flattenedData = [item for sublist in scrapedData for item in sublist]
#   tableData = list(map(list, zip(*flattenedData)))
  
#   fig = go.Figure(data=[go.Table(
#     header=dict(values=header,
#                 fill_color='paleturquoise',
#                 align='left'),
#     cells=dict(values=tableData,
#                fill_color='lavender',
#                align='left'))
# 	])
  
#   pio.show(fig)
#   #pio.write_html(fig, 'test.html')

# if __name__ == "__main__":
#   makeTable(data)

