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
        self.createFlag=0
        self.version = "0.1.0" #Major, Minor, Patch
        self.root = Tk()
        self.root.title(string="TKsXML version {}".format(self.version))
        #Inventory
        self.inv = read_inventory(inv_path)
        self.stations = [sta.code for sta in self.inv[0]]
        self.locations = []
        self._inv_to_tree(self.inv)
        #Control Board
        self._control_board()
        self.root.mainloop()
        
    def _inv_to_tree(self, inv):
        """
        Display the current StationXML Inventory.
        """
        tree = Treeview(self.root, columns=("startdate", "enddate"), height=40)
        tree.heading("startdate", text="Start Date [UTC]")
        tree.heading("enddate", text="End Date [UTC]")
        for network in inv:
            tree.insert("", 0, network.code, text=network.code)
            tree.set(network.code, "startdate", network.start_date)
            tree.set(network.code, "enddate", network.end_date)
            tree.item(network.code, open=True)
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
        
        tree.grid(row=1, column=1)
        return tree

    
    def _control_board(self):
        """
        Use a Ttk Frame widget to hold the Button commands for:
            -   adding a station-location-channel to the current network
            -   modifying an existing station-location-channel of the existing
                network
            -   deleting an existing station-location-channel of the existing 
                network
        """
                
        def _create():
            """
            When the create button is pressed:
            """
            self.bCreate.state(["disabled"])
            self.bModify.state(["disabled"])
            self.bDelete.state(["disabled"])
            self.CodeFrame = Frame(master=self.control_board)
            self.ResponseFrame = Frame(master=self.root)
            self.MetadataFrame = Frame(master=self.root)
            
            def _cancel():
                for child in self.CodeFrame.winfo_children():
                    child.destroy()
                for child in self.CodeFrameChild.winfo_children():
                    child.destroy()
                for child in self.ResponseFrame.winfo_children():
                    child.destroy()
                for child in self.MetadataFrame.winfo_children():
                    child.destroy()
                self.bCreate.state(["!disabled"])
                self.bModify.state(["!disabled"])
                self.bDelete.state(["!disabled"])
                
            def _submit():
                """
                Append the proposed metadata changes to the StationXML Inventory.
                """
                _cancel()
                pass
                
            def generate_other_menus():
                if self.createFlag == 1:
                    self.CodeFrameChild.destroy()
                    
                self.CodeFrameChild = Frame(master=self.control_board )    
                station = stationVar.get()
                if station=="Other":
                    #location
                    locationEntryVar = StringVar(self.CodeFrameChild)
                    locationEntry = Entry(self.CodeFrameChild, textvariable=locationEntryVar)
                    locationEntry.grid(row=2, column=1)
                    #channel
                    channelEntryVar = StringVar(self.CodeFrameChild)
                    channelEntry = Entry(self.CodeFrameChild, textvariable=channelEntryVar)
                    channelEntry.grid(row=4, column=1)
                    #startdate
                    startdateEntryVar = StringVar(self.CodeFrameChild)
                    startdateEntry = Entry(self.CodeFrameChild, textvariable=startdateEntryVar)
                    startdateEntry.grid(row=6, column=1)
                    #enddate
                    enddateEntryVar = StringVar(self.CodeFrameChild)
                    enddateEntry = Entry(self.CodeFrameChild, textvariable=enddateEntryVar)
                    enddateEntry.grid(row=8, column=1)
                    
                    #Other labels
                    Label(master=self.CodeFrameChild, text="Location\n(specify):", justify="center").grid(row=1, column=1, pady=10)
                    Label(master=self.CodeFrameChild, text="Channel\n(specify):", justify="center").grid(row=3, column=1, pady=10)
                    Label(master=self.CodeFrameChild, text="Start Date\n(specify)\n[YYYY-mm-ddTHH:MM:SS]:", justify="center").grid(row=5, column=1, pady=10)
                    Label(master=self.CodeFrameChild, text="End Date\n(specify)\n[YYYY-mm-ddTHH:MM:SS]:", justify="center").grid(row=7, column=1, pady=10)
                
                else:
                    location_codes = list(set([cha.location_code for sta in self.inv.select(station=station)[0] for cha in sta]))
                    channel_codes = list(set([cha.code for sta in self.inv.select(station=station)[0] for cha in sta]))
                    start_times = list(set([str(cha.start_date) for sta in self.inv.select(station=station)[0] for cha in sta]))
                    end_times = list(set([str(cha.end_date) for sta in self.inv.select(station=station)[0] for cha in sta]))
                    
                    locationVar = StringVar(self.CodeFrameChild)
                    locationVar.set(location_codes[0])
                    channelVar = StringVar(self.CodeFrameChild)
                    channelVar.set(channel_codes[0])
                    startdateVar = StringVar(self.CodeFrameChild)
                    startdateVar.set(start_times[0])
                    enddateVar = StringVar(self.CodeFrameChild)
                    enddateVar.set(end_times[0])
                    
                    #location
                    locationMenu = OptionMenu(self.CodeFrameChild, locationVar, location_codes[0], *list(location_codes + ["Other"]))
                    locationMenu.grid(row=2, column=1)
                    locationEntryVar = StringVar(self.CodeFrameChild)
                    locationEntry = Entry(self.CodeFrameChild, textvariable=locationEntryVar)
                    locationEntry.grid(row=3, column=1)
                    #channel
                    channelMenu = OptionMenu(self.CodeFrameChild, channelVar, channel_codes[0], *list(channel_codes + ["Other"]))
                    channelMenu.grid(row=5, column=1)
                    channelEntryVar = StringVar(self.CodeFrameChild)
                    channelEntry = Entry(self.CodeFrameChild, textvariable=channelEntryVar)
                    channelEntry.grid(row=6, column=1)
                    #startdate
                    startdateMenu = OptionMenu(self.CodeFrameChild, startdateVar, start_times[0], *list(start_times + ["Other"]))
                    startdateMenu.grid(row=8, column=1)
                    startdateEntryVar = StringVar(self.CodeFrameChild)
                    startdateEntry = Entry(self.CodeFrameChild, textvariable=startdateEntryVar)
                    startdateEntry.grid(row=9, column=1)
                    #enddate
                    enddateMenu = OptionMenu(self.CodeFrameChild, enddateVar, end_times[0], *list(end_times + ["Other"]))
                    enddateMenu.grid(row=11, column=1)
                    enddateEntryVar = StringVar(self.CodeFrameChild)
                    enddateEntry = Entry(self.CodeFrameChild, textvariable=enddateEntryVar)
                    enddateEntry.grid(row=12, column=1)
                    
                    #Other labels
                    Label(master=self.CodeFrameChild, text="Location\n(if Other, specify):", justify="center").grid(row=1, column=1, pady=10)
                    Label(master=self.CodeFrameChild, text="Channel\n(if Other, specify):", justify="center").grid(row=4, column=1, pady=10)
                    Label(master=self.CodeFrameChild, text="Start Date\n(if Other, specify)\n[YYYY-mm-ddTHH:MM:SS]:", justify="center").grid(row=7, column=1, pady=10)
                    Label(master=self.CodeFrameChild, text="End Date\n(if Other, specify)\n[YYYY-mm-ddTHH:MM:SS]:", justify="center").grid(row=10, column=1, pady=10)
                
                
                Button(self.CodeFrameChild, text="Cancel", command = lambda: _cancel()).grid(row=13, column=1, pady=10)
                Button(self.CodeFrameChild, text="Append to StationXML Inventory", command = lambda: _submit()).grid(row=14, column=1)
                self.CodeFrameChild.grid(row=6, column=1)
                self.createFlag=1
                return
            
            #Station Codes
            Label(master=self.CodeFrame, text="Station\n(if Other, specify):", justify="center").grid(row=1, column=1, pady=10)
            stationVar = StringVar(self.CodeFrame)
            stationVar.set(self.stations[0])
            generate_other_menus()
            stationMenu = OptionMenu(self.CodeFrame, stationVar, self.stations[0], *list(self.stations + ["Other"]), command=lambda _:generate_other_menus())
            stationMenu.grid(row=2, column=1)
            stationEntryVar = StringVar(self.CodeFrame)
            stationEntry = Entry(self.CodeFrame, textvariable=stationEntryVar)
            stationEntry.grid(row=3, column=1)
            
            self.CodeFrame.grid(row=5, column=1)
            
            #Instrument Response
            Label(master=self.ResponseFrame, text="Instrumental Response and Calibrations", justify="center", font=("Helvetica", 12, "bold")).grid(row=1, column=1, padx=25, pady=10)

            self.ResponseFrame.grid(row=1, column=3)
            
            #Channel metadata
            Label(master=self.MetadataFrame , text="Channel Metadata", justify="center", font=("Helvetica", 12, "bold")).grid(row=1, column=1, padx=50, pady=10)

            self.MetadataFrame.grid(row=1, column=4)
            return 
        
        def _modify():
            pass
        
        def _delete():
            pass
        
        self.control_board = Frame(master=self.root, height=30, width=30)
        Label(master=self.control_board , text="Control Board", justify="center", font=("Helvetica", 12, "bold")).grid(row=1, column=1, padx=150, pady=10)
        self.bCreate = Button(master=self.control_board , text="Add", command=_create)
        self.bCreate.grid(row=2, column=1)
        self.bModify = Button(master=self.control_board , text="Modify", command=_modify)
        self.bModify.grid(row=3, column=1)
        self.bDelete = Button(master=self.control_board , text="Delete", command=_delete)
        self.bDelete.grid(row=4, column=1)
        self.control_board.grid(row=1, column=2)
        
if __name__=="__main__":
    NV = TKsXML("NV.xml")