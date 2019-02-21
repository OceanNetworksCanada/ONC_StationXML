# -*- coding: utf-8 -*-
"""
Created on Wed Aug 15 14:37:00 2018

@author: jfarrugia
"""

from obspy.core.inventory import response as rp
from obspy import Inventory
from obspy.core.inventory import Network, Station, Channel, Site
import obspy
from datetime import datetime

def create_response(inputsamplerate=20.0, scaling_factor=1.4e-8, units='M/S**2'):
    _response = rp.Response(
            instrument_sensitivity=rp.InstrumentSensitivity(1.0/scaling_factor, 
                                                            1.0, 
                                                            input_units=units, 
                                                            output_units='COUNTS'
                                                            ), 
            response_stages=[
                    rp.CoefficientsTypeResponseStage(
                            stage_sequence_number=1, 
                            stage_gain=1.0/scaling_factor, 
                            stage_gain_frequency=1.0, 
                            input_units=units, 
                            output_units='COUNTS', 
                            cf_transfer_function_type='DIGITAL', 
                            numerator=rp.FloatWithUncertaintiesAndUnit(1), 
                            denominator=rp.FloatWithUncertaintiesAndUnit(1), 
                            decimation_input_sample_rate=inputsamplerate, 
                            decimation_delay=0.0, 
                            decimation_factor=1, 
                            decimation_correction=0.0, 
                            decimation_offset=0
                            )
                    ]
            )
    if units=='M/S**2':
        _response.recalculate_overall_sensitivity()
#        _response.plot(1e-3, output='ACC', label='')
        
    return _response

def create_inv(network_code, station_code, location_code, channel_code, isr, sf, u):
    writethisinv = Inventory(
            networks = [Network(code=network_code,
                                start_date=obspy.UTCDateTime('2007-01-01'),
                                stations=[Station(code=station_code,
                                                  latitude=1, 
                                                  longitude=2, 
                                                  elevation=3,
                                                  creation_date=obspy.UTCDateTime('2007-01-01'),
                                                  site=Site(name='site'),
                                                  channels=[Channel(code=channel_code, 
                                                                    location_code = location_code,
                                                                    start_date=obspy.UTCDateTime('2007-01-01'),
                                                                    latitude = 1, 
                                                                    longitude = 2, 
                                                                    elevation = 3, 
                                                                    depth = 4, 
                                                                    response=create_response(inputsamplerate=isr, scaling_factor=sf, units=u)
                                                                    )])])], 
            source='Joseph Farrugia, Ocean Networks Canada', # The source should be the id whoever create the file.
            created = obspy.UTCDateTime(datetime.today())
                )
    return writethisinv

#%%write inventories
MEMS_KEMF = ['HNE', 'HNN', 'HNZ', 'MNE', 'MNN', 'MNZ', 'HNE', 'HNN', 'HNZ', 'MNE', 'MNN', 'MNZ', 'HNE', 'HNN', 'HNZ', 'MNE', 'MNN', 'MNZ']
MEMS_KEMF_loc = ['Z1', 'Z1', 'Z1', 'Z1', 'Z1', 'Z1', 'Z2', 'Z2', 'Z2', 'Z2', 'Z2', 'Z2', 'Z3', 'Z3', 'Z3', 'Z3', 'Z3', 'Z3']
MEMS_KEMF_SR = ['100', '100', '100', '5', '5', '5', '100', '100', '100', '5', '5', '5', '100', '100', '100', '5', '5', '5']

for c, l, sr in zip(MEMS_KEMF, MEMS_KEMF_loc, MEMS_KEMF_SR):
    inv = create_inv(network_code='NV', station_code='KEMF', location_code='l', channel_code='c', isr=sr, sf=0.00059841, u='M/S**2')
    inv.write(r'\\onc-fileserver\redirect4\jfarrugia\Documents\GitHub\ONC_StationXML\{}.xml'.format("KEMF"+"."+l+"."+c+"_"+sr), format='stationxml')