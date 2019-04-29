# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 07:58:37 2018

@author: jfarrugia

This writes directly to the ONC_StationXML Git Desktop respository folder in Documents.
Changes need to be commited, and pushed to the online respository.

Merges the NV_ONC.xml file written by NV_StationXML.py with the original StationXML
    downloaded from IRIS for the NEPTUNE seismic network (NV). In this way, 
    channels for a station not accounted for in _metadata.yaml are not missed when 
    the updated metadata StationXML file is returned to IRIS.

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
                    for rs in channel.response.response_stages:
                        #correct for COUNTS units name
                        if rs.output_units=="COUNTS" or rs.output_units=="COUNT":
                            rs.output_units = rs.output_units.lower()
                            rs.output_units_description = "Digital Counts"
                        if rs.input_units=="COUNTS" or rs.output_units=="COUNT":
                            rs.input_units = rs.input_units.lower()
                            rs.input_units_description = "Digital Counts"
                        #correct for Celsius (C) units name
                        if rs.output_units=="C":
                            rs.output_units = "CELSIUS" #will be degC in the future
                        if rs.input_units=="C":
                            rs.input_units = "CELSIUS" #will be degC in the future
                    
                    if channel.response.instrument_sensitivity.output_units=="COUNTS" or channel.response.instrument_sensitivity.output_units == "COUNT":
                        channel.response.instrument_sensitivity.output_units = channel.response.instrument_sensitivity.output_units.lower()
                    if channel.response.instrument_sensitivity.input_units=="COUNTS" or channel.response.instrument_sensitivity.input_units == "COUNT":
                        channel.response.instrument_sensitivity.input_units = channel.response.instrument_sensitivity.input_units.lower()
                    if channel.calibration_units=="COUNTS" or channel.calibration_units=="COUNT": 
                        channel.calibration_units = channel.calibration_units.lower()
                    new_inv[0][k].channels.append(channel)
            
print(new_inv)

#%% Write the Inventory to StationXML
print("Writing file.")
new_inv.write(r"MetadataDrop\NV.xml", format="stationxml", validate=True)

print("\n\nStationXML is valid? {}.".format(validate_stationxml('NV_StationXML.xml')[0]))
if validate_stationxml('NV_StationXML.xml')[1] == ():
    print("\t - No errors were found.")
else:
    print("Errors found: {}".format(validate_stationxml('NV_StationXML.xml')[1]))