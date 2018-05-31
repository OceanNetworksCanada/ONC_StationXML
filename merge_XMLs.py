# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 07:58:37 2018

@author: jfarrugia

This writes directly to the ONC_StationXML Git Desktop respository folder in Documents.
Changes need to be commited, and pushed to the online respository. 

"""

from obspy import read_inventory
new_inv = read_inventory("NV_ONC.xml") #latest xml
old_inv = read_inventory("NV_stationxml-iris-download.xml") #latest xml with IRIS
from obspy.io.stationxml.core import validate_stationxml

#list of stations:
for i, stn in enumerate(old_inv[0].get_contents()['stations']):
    stn = stn.split()[0].split('.')[1]
    oldcha = old_inv.select(station=stn)[0][0]
    newcha = new_inv.select(station=stn)[0][0]
    for i, channel in enumerate(oldcha):
        if oldcha.get_contents()['channels'][i] not in newcha.get_contents()['channels']:
            newstas = new_inv[0].get_contents()['stations']
            for k in range(len(newstas)):
                if newstas[k].split()[0].split('.')[1] == stn:
                    new_inv[0][k].channels.append(channel)
            
print(new_inv)

#%% Write the Inventory to StationXML
print("Writing file.")
new_inv.write("NV_StationXML.xml", format="stationxml", validate=True)

print("\n\nStationXML is valid? {}.".format(validate_stationxml('NV_StationXML.xml')[0]))
if validate_stationxml('NV_StationXML.xml')[1] == ():
    print("\t - No errors were found.")
else:
    print("Errors found: {}".format(validate_stationxml('NV_StationXML.xml')[1]))