# -*- coding: utf-8 -*-
"""
Created on Tue Jun  4 12:00:04 2019

@author: jfarrugia
"""

from obspy import read_inventory

I = read_inventory(r"_dataloggerRESP\RESP.XX.NS000..BHZ.UNITY.DC.1.txt")