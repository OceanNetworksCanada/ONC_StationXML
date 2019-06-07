# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 10:36:26 2019

@author: Joseph Farrugia

"""

from tkinter import *
from tkinter.ttk import *
from obspy import read_inventory
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class TKsXML(object):
    
    def __init__(self, inv_path):
        self.root = Tk()
        self.inv = read_inventory(inv_path)
        tree = self._inv_to_tree(self.inv)
        tree.grid(row=1, column=1)
        self.root.mainloop()
        
    def _inv_to_tree(self, inv):
        tree = Treeview(self.root, columns=("startdate", "enddate"), height=30)
        tree.heading("startdate", text="Start Date [UTC]")
        tree.heading("enddate", text="End Date [UTC]")
        for network in inv:
            tree.insert("", 0, network.code, text=network.code)
            tree.set(network.code, "startdate", network.start_date)
            tree.set(network.code, "enddate", network.end_date)
            for station in inv.select(network=network.code)[0]:
                tree.insert(network.code, "end", station.code, text=station.code)
                tree.set(station.code, "startdate", station.start_date)
                tree.set(station.code, "enddate", station.end_date)
                for c, channel in enumerate(inv.select(station=station.code)[0][0]):
                    tree.insert(station.code, "end", "{}.{}.{}_{}".format(station.code, channel.location_code, channel.code, c), text="{}.{}.{}".format(
                            station.code, 
                            channel.location_code,
                            channel.code))
                    tree.set("{}.{}.{}_{}".format(station.code, channel.location_code, channel.code, c), "startdate", channel.start_date)
                    tree.set("{}.{}.{}_{}".format(station.code, channel.location_code, channel.code, c), "enddate", channel.end_date)

        return tree
        
        
if __name__=="__main__":
    NV = TKsXML("NV.xml")