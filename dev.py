# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 10:58:18 2019

@author: jfarrugia
"""

from tkinter import *
from tkinter import filedialog, messagebox, _setit
from tkinter.ttk import *
from obspy import read_inventory, UTCDateTime
from geopy import distance
from obspy.core.inventory import Inventory, Network, Station, Channel, Site, Equipment
from obspy.core.inventory.util import ExternalReference
import numpy as np
import os

class App(object):
    """
    An application for maintaining seismic network metadata.
    """
    
    def __init__(self, inventory_path):
        self.root = Tk()
        
        self.version = "0.1.0" #Major, Minor, Patch
        self.root.title(string="Inventory Manager v{}".format(self.version))
        
        messagebox.showinfo(title="Get Started", message="Open an Inventory file (e.g., StationXML) to get started.")
        
        #Create File Menu for Application
        self.create_menubar()
        self.datalogger_response_labels=[]
        self.datalogger_response_entries=[]
        self.sensor_response_labels=[]
        self.sensor_response_entries=[]
        
        self.inventory_name = ""
        
        os.makedirs("_dataloggerRESP", exist_ok=True)
        os.makedirs("_sensorRESP", exist_ok=True)
        os.makedirs("Inventories", exist_ok=True)
    
    def create_menubar(self):
        """
        Create the main and sub menubar for the application.
        """
        
        self.menubar = Menu(self.root)
        
        self.filemenu = Menu(self.menubar, tearoff=0)
        self.filemenu.add_command(label="Open Inventory...", command=self.open_inventory)
        self.filemenu.add_command(label="Save Inventory", command=self.save_inventory)
        self.filemenu.add_command(label="Save As Inventory...", command=self.saveas_inventory)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Selective Inventory Save...", command=self.selective_saveas_inventory)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Clear Inventory", command=self.clear_inventory)
        self.filemenu.add_separator()
        self.filemenu.add_command(label="Exit", command=self.root.destroy)
        
        self.editmenu = Menu(self.menubar, tearoff=0)
        self.editmenu.add_command(label="Append to Inventory...", command=self.append_to_inventory)
        self.editmenu.add_command(label="Modify Existing Channel...", command=self.edit_inventory_entry)
        self.editmenu.add_command(label="Remove Existing Channel...", command=self.remove_inventory_entry)
        
        self.menubar.add_cascade(label="File", menu=self.filemenu)
        self.menubar.add_cascade(label="Edit", menu=self.editmenu)
        self.menubar.entryconfig("Edit", state="disabled")
        self.root.config(menu=self.menubar)
        
        
    
    def open_inventory(self, I=None):
        """
        Function when "Open Inventory" is pressed:
        """
        if I==None:
            inventory_path = filedialog.askopenfilename(
                    initialdir = "/", 
                    title="Select Metadata Inventory", 
                    filetypes=(("xml files", ".xml"), ("all files", "*.*")))
            if inventory_path == None:
                return
            else:
                self.inv = read_inventory(inventory_path)
        
        else:
            self.inv = I
        
        self.network_tree = Treeview(self.root, columns=("startdate", "enddate"), height=40)
        self.network_tree.heading("startdate", text="Start Date [UTC]")
        self.network_tree.heading("enddate", text="End Date [UTC]")
        
        self.inventory_name = inventory_path.split("/")[-1]
        
        self.stations = [sta.code for sta in self.inv[0]]
        self.locations = []
        
        for network in self.inv:
            self.network_tree.insert("", 0, network.code, text=network.code)
            self.network_tree.set(network.code, "startdate", network.start_date)
            self.network_tree.set(network.code, "enddate", network.end_date)
            self.network_tree.item(network.code, open=True)
            for station in self.inv.select(network=network.code)[0]:
                self.network_tree.insert(network.code, "end", station.code, text=station.code)
                self.network_tree.set(station.code, "startdate", station.start_date)
                self.network_tree.set(station.code, "enddate", station.end_date)
                for c, channel in enumerate(self.inv.select(station=station.code)[0][0]):
                    self.network_tree.insert(station.code, "end", "{}.{}.{}_{}".format(station.code, channel.location_code, channel.code, c), text="{}.{}.{}".format(
                            station.code, 
                            channel.location_code,
                            channel.code))
                    self.network_tree.set("{}.{}.{}_{}".format(station.code, channel.location_code, channel.code, c), "startdate", channel.start_date)
                    self.network_tree.set("{}.{}.{}_{}".format(station.code, channel.location_code, channel.code, c), "enddate", channel.end_date)
        
        self.network_tree.pack()
        self.menubar.entryconfig("Edit", state="normal")
        return
    
    def save_inventory(self):
        
            
        if self.inventory_name == "":
            return
        else:
            self.inv.write("Inventories/" + self.inventory_name, format="STATIONXML")
            messagebox.showinfo(title="Save", 
                                    message="Success! File saved: {}".format(self.inventory_name))
        
    def saveas_inventory(self):
        if self.inventory_name == "":
            return
        else:
            filename = filedialog.asksaveasfilename(
                    initialdir = "Inventories", 
                    title="Select response file", 
                    initialfile=self.inventory_name,
                    filetypes=(("xml files", ".xml"), ("all files", "*.*")))
            if filename == None:
                return
            else:
                self.inventory_name = filename.split("/")[-1]
                self.inv.write("Inventories/" + self.inventory_name, format="STATIONXML")
                messagebox.showinfo(title="Save As", 
                                    message="Success! File saved as: {}".format(self.inventory_name))
    
    def selective_saveas_inventory(self):
        
        def selective_saveas(station, location, channel, starttime, endtime):
            
            if endtime != "*":
                try:
                    endtime = UTCDateTime(endtime)
                except:
                    messagebox.showinfo(title="Invalid End Time", 
                                        message="Invalid 'endtime' provided. Must at least provide the year-month-day (e.g., 2010-01-01).")
                    return
                
            if starttime != "*":
                try:
                    starttime = UTCDateTime(starttime)
                except:
                    messagebox.showinfo(title="Invalid Start Time", 
                                        message="Invalid 'starttime' provided. Must at least provide the year-month-day (e.g., 2010-01-01).")
                    return
                    
            sel_inv = self.inv.select(station=station, 
                                      location=location, 
                                      channel=channel, 
                                      starttime=starttime, 
                                      endtime=endtime)
            
            if sel_inv.networks == []: #empty
                messagebox.showinfo(title="Empty Inventory", 
                                    message="Parameter selection resulted in empty Inventory! \n\n{}\nStation: {}\nLocation: {}\nChannel: {}\nStart Time: {}\nEnd Time: {}\n\nFile NOT saved.".format(self.inventory_name, station, location, channel, starttime, endtime))
                return
            else:
                filename = filedialog.asksaveasfilename(
                        initialdir = "Inventories", 
                        title="Select response file", 
                        initialfile=self.inventory_name,
                        filetypes=(("xml files", ".xml"), ("all files", "*.*")))
                if filename == None:
                    return
                else:
                    self.inventory_name = filename.split("/")[-1]
                    sel_inv.write("Inventories/" + self.inventory_name, format="STATIONXML")
                    messagebox.showinfo(title="Save As", 
                                        message="Success! File saved as: {}\nStation: {}\nLocation: {}\nChannel: {}\nStart Time: {}\nEnd Time: {}".format(self.inventory_name, station, location, channel, starttime, endtime))
            
        if self.inventory_name == "":
            return
        
        else:
            popup = Toplevel()
            popup.grab_set()
            popup.title("Append to Inventory")
            
            lf = LabelFrame(master=popup, text="Parameter Selection")
            lf.pack()
            
            Label(lf, text="Station: ", justify="right").grid(row=1, column=1, sticky="e", padx=2, pady=5)
            stationEntry = Entry(master=lf)
            stationEntry.grid(row=1, column=2, sticky="w")
            stationEntry.insert(END, "*")
            
            Label(lf, text="Location: ", justify="right").grid(row=2, column=1, sticky="e", padx=2, pady=5)
            locationEntry = Entry(master=lf)
            locationEntry.grid(row=2, column=2, sticky="w")
            locationEntry.insert(END, "*")
            
            Label(lf, text="Channel: ", justify="right").grid(row=3, column=1, sticky="e", padx=2, pady=5)
            channelEntry = Entry(master=lf)
            channelEntry.grid(row=3, column=2, sticky="w")
            channelEntry.insert(END, "*")
            
            Label(lf, text="Start Time: ", justify="right").grid(row=4, column=1, sticky="e", padx=2, pady=5)
            starttimeEntry = Entry(master=lf)
            starttimeEntry.grid(row=4, column=2, sticky="w")
            starttimeEntry.insert(END, "*")
        
            Label(lf, text="End Time: ", justify="right").grid(row=5, column=1, sticky="e", padx=2, pady=5)
            endtimeEntry = Entry(master=lf)
            endtimeEntry.grid(row=5, column=2, sticky="w")
            endtimeEntry.insert(END, "*")
            
            b = Button(master=lf, text="Save As", command=lambda:selective_saveas(station=stationEntry.get(), 
                                                                              location=locationEntry.get(), 
                                                                              channel=channelEntry.get(), 
                                                                              starttime=starttimeEntry.get(), 
                                                                              endtime=endtimeEntry.get()))
            b.grid(row=6, column=1, sticky="w", padx=2, pady=5)
    
    def clear_inventory(self):
        if self.inventory_name == "":
            return
        self.network_tree.destroy()
        self.menubar.entryconfig("Edit", state="disabled")
        self.inventory_name = ""
    
    def append_to_inventory(self):

        def openResponseFile(e, _type="sensor"):
            e.delete(0, "end")
            
            try:
                if _type == "datalogger":
                    for E in self.datalogger_response_entries:
                        E.destroy()
                        self.datalogger_response_entries = []
                    for L in self.datalogger_response_labels:
                        L.destroy()
                        self.datalogger_response_labels = []
                else:
                    for E in self.sensor_response_entries:
                        E.destroy()
                        self.sensor_response_entries = []
                    for L in self.sensor_response_labels:
                        L.destroy()
                        self.sensor_response_labels = []
            except:
                pass
            
            filename = filedialog.askopenfilename(
                    initialdir = "_{}RESP".format(_type), 
                    title="Select response file", 
                    filetypes=(
                            ("xml files", ".xml"), 
                            ("txt files", ".txt"), 
                            ("sac files", ".sac"), 
                            ("all files", "*.*")))
            e.insert(END, filename)
            
            resp = read_inventory(filename)[0][0][0].response
            #starting at row 6
            if _type=="datalogger":
                move_col = 2
            else:
                move_col = 0
            
            num_response_stages = len(resp.response_stages)
            if num_response_stages > 5:
                num_response_stages = "5 of {}".format(len(resp.response_stages))
                
            l = Label(master=responseFrame,
                  justify="left",
                  text="{} response file:\n\tFrom {} to {}\n\tOverall sensitivity: {} at {} Hz\n\t{} stages:".format(
                      _type.upper(),
                      resp.instrument_sensitivity.input_units, 
                      resp.instrument_sensitivity.output_units, 
                      resp.instrument_sensitivity.value, 
                      resp.instrument_sensitivity.frequency, 
                      num_response_stages
                      ))
            l.grid(column=1 + move_col, row=6, columnspan=2, sticky="w")
            
            if _type == "sensor":
                self.sensor_response_labels.append(l)
            else:
                self.datalogger_response_labels.append(l)
                
            for i, stage in enumerate(resp.response_stages):
                if i < 5:
                    l = Label(master=responseFrame, 
                          justify="left",
                          text="\t\tStage {}: {} from {} to {}, gain: ".format(
                          i+1, 
                          stage.__str__().split("Response type: ")[1].split(", Stage")[0], 
                          stage.input_units, 
                          stage.output_units
                          ))
                    l.grid(row=7+i, column=1 + move_col, sticky="w")
                    
                    e = Entry(master=responseFrame, 
                                  width=11)
                    e.insert(END, stage.stage_gain)
                    e.grid(row=7+i, column=2 + move_col, sticky="w")
                    
                    if _type=="sensor":
                        self.sensor_response_entries.append(e)
                        self.sensor_response_labels.append(l)
                    else:
                        self.datalogger_response_entries.append(e)
                        self.datalogger_response_labels.append(l)
                
            return
        
        def disableDataloggerDialogBox(e, b):
            e.delete(0, "end")
            try:
                for E in self.datalogger_response_entries:
                    E.destroy()
                    self.datalogger_response_entries = []
                for L in self.datalogger_response_labels:
                    L.destroy()
                    self.datalogger_response_labels = []
            except:
                pass
            
            if "disabled" in e.state():
                e.state(["!disabled"])
                b.state(["!disabled"])
            else:
                e.state(["disabled"])
                b.state(["disabled"])
                
        def noCalibration(_type="sensor"):
            if _type == "sensor":
                try:
                    for E in self.sensor_response_entries:
                        if "disabled" in E.state():
                            E.state(["!disabled"])
                        else:
                            E.state(["disabled"])
                except:
                    pass
            else:
                try: 
                    for E in self.datalogger_response_entries:
                        if "disabled" in E.state():
                            E.state(["!disabled"])
                        else:
                            E.state(["disabled"])
                except:
                    pass
        
        def collate_and_parse_to_inventory(
                        network, station, location, channel, epochStartTime, epochEndTime,
                        siteName, latitude, longitude, elevation, depth, dip, azimuth, sampling_rate,
                        sensor_response, datalogger_response, instrumentType, sensorDescription, 
                        dataloggerDescription, sensorSerial, dataloggerSerial, deviceID, website):
            """
            Collate and parse the provided information into the Inventory.
            """
            def make_channel_response(_sensor_resp, _dl_resp, sensorSerial):
                _dl_resp.response_stages.pop(0)
                _dl_resp.response_stages.insert(0, _sensor_resp.response_stages[0])
                _dl_resp.instrument_sensitivity.input_units = _sensor_resp.instrument_sensitivity.input_units
                _dl_resp.instrument_sensitivity.input_units_description = _sensor_resp.instrument_sensitivity.input_units_description
                _response = _dl_resp
                
                if _response.instrument_sensitivity.output_units=="COUNTS" or _response.instrument_sensitivity.output_units == "COUNT":
                    _response.instrument_sensitivity.output_units = _response.instrument_sensitivity.output_units.lower()
                if _response.instrument_sensitivity.input_units=="COUNTS" or _response.instrument_sensitivity.input_units == "COUNT":
                    _response.instrument_sensitivity.input_units = _response.instrument_sensitivity.input_units.lower()
                    
                for stage in _response.response_stages:
                    #correct for COUNTS units name
                    if stage.output_units=="COUNTS" or stage.output_units=="COUNT":
                        stage.output_units = stage.output_units.lower()
                        stage.output_units_description = "Digital Counts"
                    if stage.input_units=="COUNTS" or stage.output_units=="COUNT":
                        stage.input_units = stage.input_units.lower()
                        stage.input_units_description = "Digital Counts"
                    #correct for Celsius (C) units name
                    if stage.output_units=="C":
                        stage.output_units = "degC"
                    if stage.input_units=="C":
                        stage.input_units = "degC"
                        
                #special condition for SA ULN 40 Vpg
                if sensorSerial.endswith("40Vpg"):
                    #ensure this change is made to an unused response stage
                    stop = 0
                    for stage in _response.response_stages:
                        if stage.stage_gain == 1 and stage.input_units == 'counts' and stage.output_units == 'counts' and stop == 0:
                            stage.stage_gain *= float(2/3)
                            stop = 1
                return _response

            # inventories can be merged using '+'
            # i.e. inv1 + inv2, or just build the station and append it
            
            inv = self.inv
            if station in self.stations: #then we just need to attach the channel to the exisiting station
                station_INDEX = np.where(np.asarray([i.code for i in inv[0].stations]) == station)[0][0]
            else: #we need to create a new station, and attach the channel
                _station = Station(
                        code = station,
                        latitude = latitude,
                        longitude = longitude, 
                        elevation = elevation, 
                        start_date = epochStartTime,
                        creation_date = None,
                        end_date = epochEndTime,
                        site = Site(name=siteName),
                        geology = None,
                        description = None
                        )
                inv[0].stations.append(_station)
                
            #calculate the response
            channel_response = make_channel_response(sensor_response, datalogger_response, sensorSerial)
            #now add the channels to the station
            _channel = Channel(
                    code = channel,
                    location_code = location,
                    latitude = latitude,
                    longitude = longitude,
                    elevation = elevation,
                    depth = depth,
                    start_date = epochStartTime,
                    end_date = epochEndTime,
                    azimuth = azimuth,
                    dip = dip,
                    types = [instrumentType],
                    sample_rate = sampling_rate,
                    equipment = Equipment(description = "{}/{}".format(sensorDescription, dataloggerDescription), serial_number = "{}/{}".format(sensorSerial, dataloggerSerial)),
                    sensor = Equipment(description = "{}/{}".format(sensorDescription, dataloggerDescription), serial_number = "{}/{}".format(sensorSerial, dataloggerSerial)),
                    response = channel_response,
                    external_references = [ExternalReference(website, 'Data Search URL.'),
                                           ExternalReference("https://data.oceannetworks.ca/DeviceListing?DeviceId=" + deviceID, 'Device URL.')]
                    )
            if _channel.calibration_units=="COUNTS" or _channel.calibration_units=="COUNT": 
                _channel.calibration_units = _channel.calibration_units.lower()
                
            if station in self.stations:
                inv[0][station_INDEX].channels.append(_channel)
            else:
                _station.channels.append(_channel)
                
            return inv
        
        def appendToInventory():
            errors_to_report = []
            warnings_to_report = []
            
            #get metadata from options menu widgets
            network = self.inv[0].code
            station = stationVar.get()
            location = locationVar.get()
            channel = channelVar.get()
            deploymentDates = deploymentDatesVar.get()
            
            #get metadata from entry widgets
            siteName = siteNameEntry.get()
            latitude = latitudeEntry.get()
            longitude = longitudeEntry.get()
            elevation = elevationEntry.get()
            sampleRate = sampleRateEntry.get()
            azimuth = azimuthEntry.get()
            depth = depthEntry.get()
            dip = dipEntry.get()
            
            #response information
            sensor_response_file = sensorResponseEntry.get()
            if sensor_response_file == "":
                errors_to_report.append("- Choose a SENSOR response file.")
                sensor_response = None
            else:
                try:
                    sensor_response = read_inventory(sensor_response_file)[0][0][0].response
                except:
                    errors_to_report.append("- Error reading in chosen SENSOR response file.")
               
            datalogger_response_file = dataloggerResponseEntry.get()    
            if datalogger_response_file == "" and "disabled" not in dataloggerResponseEntry.state():
                errors_to_report.append("- Choose a DATALOGGER response file.")
            elif "disabled" in dataloggerResponseEntry.state():
                datalogger_response = sensor_response
            else:
                try:
                    datalogger_response = read_inventory(datalogger_response_file)[0][0][0].response
                except:
                    errors_to_report.append("- Error reading in chosen DATALOGGER response file.")
                
            if station=="Other": #get metadata from entry widgets
                station = stationEntry.get().strip().upper()
                
            if location == "Other" or "disabled" in locationMenu.state():
                location = locationEntry.get().strip().upper()
                
            if channel == "Other" or "disabled" in channelMenu.state():
                channel = channelEntry.get().strip().upper()
            
            if deploymentDates == "Other" or "disabled" in deploymentDatesMenu.state():
                try:
                    deployment_times_from_entry = deploymentDatesEntry.get()
                    epochStartTime = UTCDateTime(deployment_times_from_entry.split("--")[0].strip())
                    epochEndTime = deployment_times_from_entry.split("--")[1].strip()
                    if epochEndTime != "" and epochEndTime != "None":
                        epochEndTime = UTCDateTime(epochEndTime)
                        if epochStartTime > epochEndTime:
                            errors_to_report.append("- The epoch start date exceeds the epoch end time")
                except:
                    errors_to_report.append("- Provide a deployment in the format 'starttime -- endtime', where starttime is 'YYYY-mm-ddTHH:MM:SS', and endtime can be expressed as 'YYYY-mm-ddTHH:MM:SS' or 'None'")
            else:
                epochStartTime = UTCDateTime(deploymentDates.split("--")[0].strip())
                epochEndTime = deploymentDates.split("--")[1].strip()
                if epochEndTime == "" or epochEndTime == "None":
                    epochEndTime = "None"
                else: 
                    epochEndTime = UTCDateTime(epochEndTime)
                   
                    
            if (len(station) > 5 or station=="") or station.isalnum() == False:
                errors_to_report.append("- Provide an alphanumeric station code no more than five characters long")
                
            if location != "" and len(location) < 2 and location.isalnum() == False:
                errors_to_report.append("- Provide an alphanumeric location code two characters in length, or leave blank")
                
            if channel == "" or channel.isalnum() == False:
                errors_to_report.append("- Provide an alphanumeric channel code two to three characters in length")                
            
            
            #check if the proposed channel conflicts with contents of inventory
            proposed_channel = "{}.{}.{}.{}".format(network, station, location, channel)
            if proposed_channel in self.inv.get_contents()["channels"]:
                #check the deployment times
                try:
                    if epochEndTime == "None":
                        _endTimeCheck = UTCDateTime("2599-01-01")
                    else:
                        _endTimeCheck = UTCDateTime(epochEndTime)
                        
                    search_channel = self.inv.select(
                            network=network,
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=epochStartTime, 
                            endtime=_endTimeCheck)[0][0][0]
                    if search_channel.end_date == None: #should modify first 
                        errors_to_report.append("- {} with Epoch Start Date: {} and Epoch End Date: {} conflicts with existing channel with encompassing deployment dates: \n\n{}\n\n- > Modify the existing channel first, then append the new channel.".format(proposed_channel, epochStartTime, epochEndTime, search_channel))
                    else:
                        errors_to_report.append("- {} with Epoch Start Date: {} and Epoch End Date: {} conflicts with existing channel with encompassing deployment dates: \n\n{}".format(proposed_channel, epochStartTime, epochEndTime, search_channel))
                except: 
                    pass
            
            #distance between stations criteria
            if station in self.stations:
                try:
                    inv_latitude = self.inv.select(station=station, location="")[0][0][0].latitude
                    inv_longitude = self.inv.select(station=station, location="")[0][0][0].longitude
                except:
                    inv_latitude = self.inv.select(station=station)[0][0][0].latitude
                    inv_longitude = self.inv.select(station=station)[0][0][0].longitude
                    
                _inv_station_location = (inv_latitude, inv_longitude)
                _app_station_location = (latitude, longitude)
                station_separation = distance.geodesic(_inv_station_location, _app_station_location).km
                
                if station_separation > 1.5: #km
                    errors_to_report.append("- Stations with the same code are separated by a distance greater than 1.5 km.")
            else:
                _inv_station_locations = list(set([(cha.latitude, cha.longitude) for sta in self.inv[0] for cha in sta]))
                _app_station_location = (latitude, longitude)
                distances = list(sorted(set([distance.geodesic(_app_station_location, l).km for l in _inv_station_locations])))
                if distances[0] <= 1.5:
                    errors_to_report.append("- Provided (latitude, longitude) pair <= 1.5 km from station with a different code.")
            
            if siteName == "":
                warnings_to_report.append("- Site name not provided")
                
            if latitude == "":
                errors_to_report.append("- latitude not provided")
                
            if longitude == "":
                errors_to_report.append("- longitude not provided")
                
            try:
                elevation = float(elevation)
                    
            except ValueError:
                errors_to_report.append("- Provide a numeric elevation [masl]")
                
            try:
                sampleRate = float(sampleRate)
                if sampleRate < 0:
                    errors_to_report.append("- Provide a numeric sample rate >= 0")
                elif sampleRate == 0:
                    warnings_to_report.append("- A sampling rate of '0' was entered. If not an error, ignore.")
                    
            except ValueError:
                errors_to_report.append("- Provide a numeric sample rate >= 0")
                
            try:
                azimuth = float(azimuth)
                if azimuth < 0 or azimuth >= 360:
                    errors_to_report.append("- Provide a numeric azimuth (0 <= azimuth < 360)")
                    
            except ValueError:
                errors_to_report.append("- Provide a numeric azimuth (0 <= azimuth < 360)")
                
            if channel.endswith("N") and (azimuth < 355 and azimuth > 5):
                warnings_to_report.append("- Azimuth of {} was given for N-component channel, which is typically within 5 degrees of True North.".format(azimuth))
            elif channel.endswith("1") and (azimuth >= 355 or azimuth <= 5):
                warnings_to_report.append("- Channel orientation code '1' but azimuth ({}) is within 5 degress of True North so should end with 'N'".format(azimuth))
            
            if channel.endswith("E") and (azimuth < 85 or azimuth > 95):
                warnings_to_report.append("- Azimuth of {} was given for E-component channel, which is typically within 5 degrees of East.".format(azimuth))
            elif channel.endswith("2") and (azimuth >= 85 or azimuth <= 95):
                warnings_to_report.append("- Channel orientation code '2' but azimuth ({}) is within 5 degress of East so should end with 'E'".format(azimuth))
            
            samplerate = sampleRate
            if channel[0]=="E" and (samplerate < 80 or samplerate > 250):
                warnings_to_report.append("- Sample rate ({}) for 'extremely short period' channel outside normal bounds (>= 80 to < 250)".format(samplerate))
            elif channel[0]=="S" and (samplerate < 10 or samplerate > 80):
                warnings_to_report.append("- Sample rate ({}) for 'short period' channel outside normal bounds (>= 10 to < 80)".format(samplerate))
            elif channel[0]=="H" and (samplerate < 80 or samplerate > 250):
                warnings_to_report.append("- Sample rate ({}) for 'high broad band' channel outside normal bounds (>= 80 to < 250)".format(samplerate))
            elif channel[0]=="B" and (samplerate < 10 or samplerate > 80):
                warnings_to_report.append("- Sample rate ({}) for 'broad band' channel outside normal bounds (>= 10 to < 80)".format(samplerate))
            elif channel[0]=="M" and (samplerate < 1 or samplerate > 10):
                warnings_to_report.append("- Sample rate ({}) for 'mid period' channel outside normal bounds (>= 1 to < 10)".format(samplerate))
            elif channel[0]=="L" and samplerate != 1:
                warnings_to_report.append("- Sample rate ({}) for 'long period' channel outside normal bounds (~ 1)".format(samplerate))
            elif channel[0]=="V" and samplerate != 0.1:
                warnings_to_report.append("- Sample rate ({}) for 'very long period' channel outside normal bounds (~ 0.1)".format(samplerate))
            elif channel[0]=="U" and samplerate != 0.01:
                warnings_to_report.append("- Sample rate ({}) for 'ultra long period' channel outside normal bounds (~ 0.01)".format(samplerate))
            elif channel[0]=="R" and (samplerate < 0.0001 or samplerate > 0.001):
                warnings_to_report.append("- Sample rate ({}) for 'extremely long period' channel outside normal bounds (>= 0.0001 to < 0.001)".format(samplerate))
            elif channel[0]=="P" and (samplerate < 0.00001 or samplerate > 0.0001):
                warnings_to_report.append("- Sample rate ({}) for '~1 day sampling' channel outside normal bounds (>= 0.00001 to < 0.0001)".format(samplerate))
            elif channel[0]=="T" and (samplerate < 0.000001 or samplerate > 0.00001):
                warnings_to_report.append("- Sample rate ({}) for '~10 day sampling' channel outside normal bounds (>= 0.000001 to < 0.00001)".format(samplerate))
            elif channel[0]=="Q" and samplerate > 0.000001:
                warnings_to_report.append("- Sample rate ({}) for '>10 day sampling' channel outside normal bounds (< 0.000001)".format(samplerate))
            elif channel[0] not in "ESHBMLVURPTQAO":
                warnings_to_report.append("- Band code in '{}' (first character) not recognized according to SEED Manual Appendix A".format(channel))
                
            try:
                depth = float(depth)
                if depth < 0:
                    errors_to_report.append("- Provide a numeric depth >= 0 [m]")
            except ValueError:
                errors_to_report.append("- Provide a numeric depth >= 0 [m]")
            
            try:
                dip = float(dip)
                if dip < -90 or dip > 90:
                    errors_to_report.append("- Provide a numeric dip (-90 <= dip <= 90)")
                    
            except ValueError:
                errors_to_report.append("- Provide a numeric dip (-90 <= dip <= 90)")
            
            if channel.endswith("Z") and (dip > -85):
                warnings_to_report.append("- Dip of {} was given for Z-component channel, which is typically within 5 degrees of vertical.".format(dip))
            elif channel.endswith("3") and (dip <= -85):
                warnings_to_report.append("- Channel orientation code is '3' but dip ({}) is within 5 degress of vertical so should end with 'Z'".format(azimuth))
            
            instrumentType = instrumentTypeEntry.get()
            if instrumentType == "":
                errors_to_report.append("- Provide an instrument type.")
            
            sensorDescription = sensorDescriptionEntry.get()
            if sensorDescription == "":
                errors_to_report.append("- Provide a SENSOR description.")
            
            dataloggerDescription = dataloggerDescriptionEntry.get()
            if dataloggerDescription == "":
                errors_to_report.append("- Provide a DATALOGGER description.")
            
            sensorSerial = sensorSerialEntry.get()
            dataloggerSerial = dataloggerSerialEntry.get()
            deviceID = deviceIDEntry.get()
            website = websiteEntry.get()
            
            #some additional response information checks
            if self.datalogger_response_entries==[] and dataloggerResponseVar.get()==0:
                errors_to_report.append("- If not loading a DATALOGGER response file, check the box to use SENSOR stage gains.")
                
            if sensorCalibrationVar.get() == 0: #not checked, may have changed gains
                for i, entry in enumerate(self.sensor_response_entries):
                    gain = entry.get()
                    try:
                        gain = float(gain)
                        sensor_response.response_stages[i].stage_gain = gain
                    except ValueError:
                        errors_to_report.append("- Enter a valid SENSOR gain value [float or int] for Stage {}.".format(i+1))
                try:
                    sensor_response.recalculate_overall_sensitivity()
                except:
                    warnings_to_report.append("- Could not recalculate overall sensitivity for SENSOR response information for channel '{}'.".format(channel))
                    
            if "disabled" not in dataloggerResponseEntry.state() and dataloggerResponseVar.get()==0:
                for i, entry in enumerate(self.datalogger_response_entries):
                    gain = entry.get()
                    try:
                        gain = float(gain)
                        datalogger_response.response_stages[i].stage_gain = gain
                    except ValueError:
                        errors_to_report.append("- Enter a valid DATALOGGER gain value [float or int] for Stage {}.".format(i+1))
                try:
                    datalogger_response.recalculate_overall_sensitivity()
                except:
                    warnings_to_report.append("- Could not recalculate overall sensitivity for DATALOGGER response information for channel '{}'.".format(channel))
                    
                
            
            #append or don't and report errors, report any warnings
            if errors_to_report == []:
                if warnings_to_report != []:
                    messagebox.showwarning(
                            title="Warnings ({})".format(len(warnings_to_report)), 
                            message="\n".join(warnings_to_report))
                
                #actually attach to the Inventory
                self.inv = collate_and_parse_to_inventory(
                        network, station, location, channel, epochStartTime, epochEndTime,
                        siteName, latitude, longitude, elevation, depth, dip, azimuth, samplerate,
                        sensor_response, datalogger_response, instrumentType, sensorDescription, 
                        dataloggerDescription, sensorSerial, dataloggerSerial, deviceID, website)
                
                self.clear_inventory()
                self.open_inventory(I=self.inv)
                
                message = """
                      Appended:
                          {}
                      """.format(self.inv.select(
                      station=station, 
                      location=location, 
                      channel=channel,
                      starttime=epochStartTime, 
                      endtime=epochEndTime))
#                      
                messagebox.showinfo(title="Appended to Inventory", 
                                    message=message)
                
            else:
                messagebox.showerror(title="Errors ({})".format(len(errors_to_report)), 
                                     message="\n".join(errors_to_report))
                
                if warnings_to_report != []:
                    messagebox.showwarning(
                            title="Warnings ({})".format(len(warnings_to_report)), 
                            message="\n".join(warnings_to_report))
            
                
        def station_code_refresh():
            station = stationVar.get()
            
            if station != "Other":
            
                location_codes = list(sorted(set([cha.location_code for sta in self.inv.select(
                        station=station)[0] for cha in sta])))
                location = location_codes[0]
                channel_codes = list(sorted(set([cha.code for sta in self.inv.select(
                        station=station, 
                        location=location)[0] for cha in sta])))
                channel = channel_codes[0]
                
                deployment_times = list(sorted(set([(str(cha.start_date), str(cha.end_date)) for sta in self.inv.select(
                        station=stationVar.get(), 
                        location=location,
                        channel=channel)[0] for cha in sta])))
    
                starttime = deployment_times[0][0]
                endtime = deployment_times[0][1]
                starttime = UTCDateTime(starttime)
                
                deployment_times_for_entry = []
                for d in deployment_times:
                    deployment_times_for_entry.append(" -- ".join(d))

                if endtime != "" and endtime != "None": 
                    endtime = UTCDateTime(endtime)
                elif endtime == "None" or endtime == "":
                    endtime = None
                    
                locationMenu['menu'].delete(0, 'end')
                for choice in location_codes + ["Other"]:
                    locationMenu['menu'].add_command(label=choice, command=_setit(
                            locationVar, 
                            choice, 
                            lambda _: location_code_refresh()
                            ))
                locationVar.set(location_codes[0])
                
                #refresh channel codes based on station choice
                channelMenu['menu'].delete(0, 'end')
                for choice in channel_codes + ["Other"]:
                    channelMenu['menu'].add_command(label=choice, command=_setit(
                            channelVar, 
                            choice,
                            lambda _: channel_code_refresh()
                            ))
                channelVar.set(channel_codes[0])
                
                #refresh the deployment dates
                deploymentDatesMenu['menu'].delete(0, 'end')
                for choice in deployment_times_for_entry + ["Other"]:
                    deploymentDatesMenu['menu'].add_command(label=choice, command=_setit(
                            deploymentDatesVar, 
                            choice,
                            lambda _: deployment_times_refresh()
                            ))
                deploymentDatesVar.set(deployment_times_for_entry[0])
                
                #Other channel parameters
                sitename = self.inv.select(
                        station=station, 
                        location=location)[0][0].site.name
                latitude = self.inv.select(
                        station=station, 
                        location=location)[0][0][0].latitude
                longitude = self.inv.select(
                        station=station, 
                        location=location)[0][0][0].longitude
                elevation = self.inv.select(
                        station=station, 
                        location=location)[0][0][0].elevation
                samplerate = self.inv.select(
                        station=station, 
                        location=location, 
                        channel=channel, 
                        starttime=starttime, 
                        endtime=endtime)[0][0][0].sample_rate
                azimuth = self.inv.select(
                        station=station, 
                        location=location, 
                        channel=channel, 
                        starttime=starttime, 
                        endtime=endtime)[0][0][0].azimuth
                dip = self.inv.select(
                        station=station, 
                        location=location, 
                        channel=channel, 
                        starttime=starttime, 
                        endtime=endtime)[0][0][0].dip
                depth = self.inv.select(
                        station=station, 
                        location=location, 
                        channel=channel, 
                        starttime=starttime, 
                        endtime=endtime)[0][0][0].depth
                
                try:
                    instrument_types = self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].types[0]
                except:
                    instrument_types = ""
                    
                try:
                    sensor_description = self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sensor.description.split("/")[0]
                except:
                    sensor_description = ""
                    
                try:
                    datalogger_description = self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sensor.description.split("/")[1]
                except:
                    datalogger_description = ""
                
                try:
                    sensor_serial_number = self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sensor.serial_number.split("/")[0]
                except:
                    sensor_serial_number = ""
                
                try:
                    datalogger_serial_number = self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sensor.serial_number.split("/")[1]
                except:
                    datalogger_serial_number = ""
                
                try:
                    external_references = self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].external_references
                    if "Data Search" in external_references[0].description:
                        website = external_references[0].uri
                    else:
                        website = external_references[1].uri
                    if "Device" in external_references[0].description:
                        device_id = external_references[0].uri.split("=")[1]
                    else:
                        device_id = external_references[1].uri.split("=")[1]
                except:
                    device_id = ""
                    website = ""
                    
                siteNameEntry.delete(0, 'end')
                latitudeEntry.delete(0, 'end')
                longitudeEntry.delete(0, 'end')
                elevationEntry.delete(0, 'end')
                sampleRateEntry.delete(0, 'end')
                azimuthEntry.delete(0, 'end')
                depthEntry.delete(0, 'end')
                dipEntry.delete(0, 'end')
                instrumentTypeEntry.delete(0, 'end')
                sensorDescriptionEntry.delete(0, 'end')
                dataloggerDescriptionEntry.delete(0, 'end')
                sensorSerialEntry.delete(0, 'end')
                dataloggerSerialEntry.delete(0, 'end')
                deviceIDEntry.delete(0, 'end')
                websiteEntry.delete(0, 'end')
                
                stationEntry.state(["disabled"])
                locationEntry.state(["disabled"])
                channelEntry.state(["disabled"])
                deploymentDatesEntry.state(["disabled"])
                locationMenu.state(["!disabled"])
                channelMenu.state(["!disabled"])
                deploymentDatesMenu.state(["!disabled"])
                
                siteNameEntry.insert(END, str(sitename))
                latitudeEntry.insert(END, str(latitude))
                longitudeEntry.insert(END, str(longitude))
                elevationEntry.insert(END, str(elevation))
                instrumentTypeEntry.insert(END, str(instrument_types))
                sensorDescriptionEntry.insert(END, str(sensor_description))
                dataloggerDescriptionEntry.insert(END, str(datalogger_description))
                sensorSerialEntry.insert(END, str(sensor_serial_number))
                dataloggerSerialEntry.insert(END, str(datalogger_serial_number))
                deviceIDEntry.insert(END, str(device_id))
                websiteEntry.insert(END, str(website))
                sampleRateEntry.insert(END, samplerate)
                azimuthEntry.insert(END, azimuth)
                depthEntry.insert(END, depth)
                dipEntry.insert(END, dip)
                
            else:
                stationEntry.state(["!disabled"])
                locationEntry.state(["!disabled"])
                channelEntry.state(["!disabled"])
                deploymentDatesEntry.state(["!disabled"])
                locationMenu.state(["disabled"])
                channelMenu.state(["disabled"])
                deploymentDatesMenu.state(["disabled"])
                
            return
        
        def location_code_refresh():
            station = stationVar.get()
            location= locationVar.get()
            
            if location=="Other":
                locationEntry.state(["!disabled"])
                
            else:
                
                #first check if other entries are set to Other
                if (channelVar.get() == "Other" or deploymentDatesVar.get() == "Other") and location != "Other":
                    locationEntry.state(["disabled"])
                    #keep disabled
                    #we can't do anything with 'other' information
                    
                
                else:
                    locationEntry.state(["disabled"])
                    
                    channel_codes = list(sorted(set([cha.code for sta in self.inv.select(
                            station=station, 
                            location=location)[0] for cha in sta])))
                    channel = channel_codes[0]
                    
                    deployment_times = list(sorted(set([(str(cha.start_date), str(cha.end_date)) for sta in self.inv.select(
                            station=station, 
                            location=location,
                            channel=channel)[0] for cha in sta])))
        
                    starttime = deployment_times[0][0]
                    endtime = deployment_times[0][1]
                    starttime = UTCDateTime(starttime)
                    
                    deployment_times_for_entry = []
                    for d in deployment_times:
                        deployment_times_for_entry.append(" -- ".join(d))
    
                    if endtime != "" and endtime != "None": 
                        endtime = UTCDateTime(endtime)
                    elif endtime == "None" or endtime == "":
                        endtime = None
                    
                    #refresh channel codes based on station choice
                    channelMenu['menu'].delete(0, 'end')
                    for choice in channel_codes + ["Other"]:
                        channelMenu['menu'].add_command(label=choice, command=_setit(
                                channelVar, 
                                choice,
                                lambda _: channel_code_refresh()
                                ))
                    channelVar.set(channel_codes[0])
                    
                    #refresh the start dates
                    deploymentDatesMenu['menu'].delete(0, 'end')
                    for choice in deployment_times_for_entry + ["Other"]:
                        deploymentDatesMenu['menu'].add_command(label=choice, command=_setit(
                                deploymentDatesVar, 
                                choice,
                                lambda _: deployment_times_refresh()
                                ))
                    deploymentDatesVar.set(deployment_times_for_entry[0])
                    
                    if location == "Other":
                        locationEntry.state(["!disabled"])
                        location = None
                    if channel == "Other":
                        channel = None
                    if starttime == "Other":
                        starttime = None
                    if endtime == "Other":
                        endtime = None
                    
                    siteNameEntry.delete(0, 'end')
                    latitudeEntry.delete(0, 'end')
                    longitudeEntry.delete(0, 'end')
                    elevationEntry.delete(0, 'end')
                    sampleRateEntry.delete(0, 'end')
                    azimuthEntry.delete(0, 'end')
                    depthEntry.delete(0, 'end')
                    dipEntry.delete(0, 'end')
                    instrumentTypeEntry.delete(0, 'end')
                    sensorDescriptionEntry.delete(0, 'end')
                    dataloggerDescriptionEntry.delete(0, 'end')
                    sensorSerialEntry.delete(0, 'end')
                    dataloggerSerialEntry.delete(0, 'end')
                    deviceIDEntry.delete(0, 'end')
                    websiteEntry.delete(0, 'end')
                    
                    try:
                        instrument_types = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].types[0]
                    except:
                        instrument_types = ""
                        
                    try:
                        sensor_description = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.description.split("/")[0]
                    except:
                        sensor_description = ""
                        
                    try:
                        datalogger_description = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.description.split("/")[1]
                    except:
                        datalogger_description = ""
                    
                    try:
                        sensor_serial_number = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.serial_number.split("/")[0]
                    except:
                        sensor_serial_number = ""
                    
                    try:
                        datalogger_serial_number = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.serial_number.split("/")[1]
                    except:
                        datalogger_serial_number = ""
                    
                    try:
                        external_references = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].external_references
                        if "Data Search" in external_references[0].description:
                            website = external_references[0].uri
                        else:
                            website = external_references[1].uri
                        if "Device" in external_references[0].description:
                            device_id = external_references[0].uri.split("=")[1]
                        else:
                            device_id = external_references[1].uri.split("=")[1]
                    except:
                        device_id = ""
                        website = ""
                        
                    siteNameEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0].site.name))
                    latitudeEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].latitude))
                    longitudeEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].longitude))
                    elevationEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].elevation))
                    sampleRateEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sample_rate))
                    azimuthEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].azimuth))
                    dipEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].dip))
                    depthEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].depth))
                    
                    instrumentTypeEntry.insert(END, str(instrument_types))
                    sensorDescriptionEntry.insert(END, str(sensor_description))
                    dataloggerDescriptionEntry.insert(END, str(datalogger_description))
                    sensorSerialEntry.insert(END, str(sensor_serial_number))
                    dataloggerSerialEntry.insert(END, str(datalogger_serial_number))
                    deviceIDEntry.insert(END, str(device_id))
                    websiteEntry.insert(END, str(website))
            return
        
        def channel_code_refresh():
            
            station = stationVar.get()
            location= locationVar.get()
            channel = channelVar.get()
            if channel=="Other":
                channelEntry.state(["!disabled"])
            
            
            else:
                
                #first check if other entries are set to Other
                if (location == "Other" or deploymentDatesVar.get() == "Other") and channelVar.get() != "Other":
                    channelEntry.state(["disabled"]) 
                    #keep disabled
                    #we can't do anything with 'other' information
                    
                else:
                    channelEntry.state(["disabled"])
                    
                    deployment_times = list(sorted(set([(str(cha.start_date), str(cha.end_date)) for sta in self.inv.select(
                            station=station, 
                            location=location,
                            channel=channel)[0] for cha in sta])))
        
                    starttime = deployment_times[0][0]
                    endtime = deployment_times[0][1]
                    starttime = UTCDateTime(starttime)
                    
                    deployment_times_for_entry = []
                    for d in deployment_times:
                        deployment_times_for_entry.append(" -- ".join(d))
    
                    if endtime != "" and endtime != "None": 
                        endtime = UTCDateTime(endtime)
                    elif endtime == "None" or endtime == "":
                        endtime = None
                    
                    #refresh the start dates
                    deploymentDatesMenu['menu'].delete(0, 'end')
                    for choice in deployment_times_for_entry + ["Other"]:
                        deploymentDatesMenu['menu'].add_command(label=choice, command=_setit(
                                deploymentDatesVar, 
                                choice,
                                lambda _: deployment_times_refresh()
                                ))
                    deploymentDatesVar.set(deployment_times_for_entry[0])
                    
                    if starttime == "Other":
                        starttime = None
                    if endtime == "Other":
                        endtime = None
                    
                    siteNameEntry.delete(0, 'end')
                    latitudeEntry.delete(0, 'end')
                    longitudeEntry.delete(0, 'end')
                    elevationEntry.delete(0, 'end')
                    sampleRateEntry.delete(0, 'end')
                    azimuthEntry.delete(0, 'end')
                    depthEntry.delete(0, 'end')
                    dipEntry.delete(0, 'end')
                    instrumentTypeEntry.delete(0, 'end')
                    sensorDescriptionEntry.delete(0, 'end')
                    dataloggerDescriptionEntry.delete(0, 'end')
                    sensorSerialEntry.delete(0, 'end')
                    dataloggerSerialEntry.delete(0, 'end')
                    deviceIDEntry.delete(0, 'end')
                    websiteEntry.delete(0, 'end')
                    
                    try:
                        instrument_types = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].types[0]
                    except:
                        instrument_types = ""
                        
                    try:
                        sensor_description = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.description.split("/")[0]
                    except:
                        sensor_description = ""
                        
                    try:
                        datalogger_description = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.description.split("/")[1]
                    except:
                        datalogger_description = ""
                    
                    try:
                        sensor_serial_number = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.serial_number.split("/")[0]
                    except:
                        sensor_serial_number = ""
                    
                    try:
                        datalogger_serial_number = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.serial_number.split("/")[1]
                    except:
                        datalogger_serial_number = ""
                    
                    try:
                        external_references = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].external_references
                        if "Data Search" in external_references[0].description:
                            website = external_references[0].uri
                        else:
                            website = external_references[1].uri
                        if "Device" in external_references[0].description:
                            device_id = external_references[0].uri.split("=")[1]
                        else:
                            device_id = external_references[1].uri.split("=")[1]
                    except:
                        device_id = ""
                        website = ""
                        
                    siteNameEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0].site.name))
                    latitudeEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].latitude))
                    longitudeEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].longitude))
                    elevationEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].elevation))
                    sampleRateEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sample_rate))
                    azimuthEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].azimuth))
                    dipEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].dip))
                    depthEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].depth))
                    
                    instrumentTypeEntry.insert(END, str(instrument_types))
                    sensorDescriptionEntry.insert(END, str(sensor_description))
                    dataloggerDescriptionEntry.insert(END, str(datalogger_description))
                    sensorSerialEntry.insert(END, str(sensor_serial_number))
                    dataloggerSerialEntry.insert(END, str(datalogger_serial_number))
                    deviceIDEntry.insert(END, str(device_id))
                    websiteEntry.insert(END, str(website))
            return
        
        def deployment_times_refresh():
            
            
            station = stationVar.get()
            location= locationVar.get()
            channel = channelVar.get()
            deployment_times = deploymentDatesVar.get()
            
            if deployment_times=="Other":
                deploymentDatesEntry.state(["!disabled"])
            
            else:
                if (locationVar.get() == "Other" or channelVar.get() == "Other") and deployment_times != "Other":
                    deploymentDatesEntry.state(["disabled"]) 
                    #keep disabled
                    #we can't do anything with 'other' information
                else:
                    deploymentDatesEntry.state(["disabled"])
                    
                    starttime = deployment_times.split("--")[0].strip()
                    endtime = deployment_times.split("--")[1].strip()
                    starttime = UTCDateTime(starttime)
    
                    if endtime != "" and endtime != "None": 
                        endtime = UTCDateTime(endtime)
                    elif  endtime == "" or endtime == "None":
                        endtime = None
                    
                    siteNameEntry.delete(0, 'end')
                    latitudeEntry.delete(0, 'end')
                    longitudeEntry.delete(0, 'end')
                    elevationEntry.delete(0, 'end')
                    sampleRateEntry.delete(0, 'end')
                    azimuthEntry.delete(0, 'end')
                    depthEntry.delete(0, 'end')
                    dipEntry.delete(0, 'end')
                    instrumentTypeEntry.delete(0, 'end')
                    sensorDescriptionEntry.delete(0, 'end')
                    dataloggerDescriptionEntry.delete(0, 'end')
                    sensorSerialEntry.delete(0, 'end')
                    dataloggerSerialEntry.delete(0, 'end')
                    deviceIDEntry.delete(0, 'end')
                    websiteEntry.delete(0, 'end')
                    
                    try:
                        instrument_types = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].types[0]
                    except:
                        instrument_types = ""
                        
                    try:
                        sensor_description = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.description.split("/")[0]
                    except:
                        sensor_description = ""
                        
                    try:
                        datalogger_description = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.description.split("/")[1]
                    except:
                        datalogger_description = ""
                    
                    try:
                        sensor_serial_number = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.serial_number.split("/")[0]
                    except:
                        sensor_serial_number = ""
                    
                    try:
                        datalogger_serial_number = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].sensor.serial_number.split("/")[1]
                    except:
                        datalogger_serial_number = ""
                    
                    try:
                        external_references = self.inv.select(
                                station=station, 
                                location=location, 
                                channel=channel, 
                                starttime=starttime, 
                                endtime=endtime)[0][0][0].external_references
                        if "Data Search" in external_references[0].description:
                            website = external_references[0].uri
                        else:
                            website = external_references[1].uri
                        if "Device" in external_references[0].description:
                            device_id = external_references[0].uri.split("=")[1]
                        else:
                            device_id = external_references[1].uri.split("=")[1]
                    except:
                        device_id = ""
                        website = ""
                        
                    siteNameEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0].site.name))
                    latitudeEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].latitude))
                    longitudeEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].longitude))
                    elevationEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location)[0][0][0].elevation))
                    sampleRateEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].sample_rate))
                    azimuthEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].azimuth))
                    dipEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].dip))
                    depthEntry.insert(END, str(self.inv.select(
                            station=station, 
                            location=location, 
                            channel=channel, 
                            starttime=starttime, 
                            endtime=endtime)[0][0][0].depth))
                    
                    instrumentTypeEntry.insert(END, str(instrument_types))
                    sensorDescriptionEntry.insert(END, str(sensor_description))
                    dataloggerDescriptionEntry.insert(END, str(datalogger_description))
                    sensorSerialEntry.insert(END, str(sensor_serial_number))
                    dataloggerSerialEntry.insert(END, str(datalogger_serial_number))
                    deviceIDEntry.insert(END, str(device_id))
                    websiteEntry.insert(END, str(website))
            return
                
        popup = Toplevel()
        popup.grab_set()
        popup.title("Append to Inventory")
        
        #Channel metadata
        channelFrame = LabelFrame(popup, text="Channel Metadata", relief=RIDGE)
        channelFrame.grid(row=1, column=1, padx=3, pady=3, sticky="w")
        
        Label(master=channelFrame, text="Station (if set to 'Other', specify)", justify="left").grid(row=1, column=1, sticky="w", padx=5, pady=2)
        
        stationVar = StringVar(channelFrame)
        locationVar = StringVar(channelFrame)
        channelVar = StringVar(channelFrame)
        deploymentDatesVar = StringVar(channelFrame)
        
        stationVar.set(self.stations[0])
        
        location_codes = list(sorted(set([cha.location_code for sta in self.inv.select(
                station=stationVar.get())[0] for cha in sta])))
        location = location_codes[0]
        channel_codes = list(sorted(set([cha.code for sta in self.inv.select(
                station=stationVar.get(), 
                location=location)[0] for cha in sta])))
        channel = channel_codes[0]
        deployment_times = list(sorted(set([(str(cha.start_date), str(cha.end_date)) for sta in self.inv.select(
                station=stationVar.get(), 
                location=location,
                channel=channel)[0] for cha in sta])))
        starttime = deployment_times[0][0]
        endtime = deployment_times[0][1]
        starttime = UTCDateTime(starttime)
        
        deployment_times_for_entry = []
        for d in deployment_times:
            deployment_times_for_entry.append(" -- ".join(d))
        
        if endtime != "None": 
            endtime = UTCDateTime(endtime)
        elif endtime == "None":
            endtime = None
                    
        stationMenu = OptionMenu(
                channelFrame, 
                stationVar, 
                self.stations[0], 
                *list(self.stations + ["Other"]), 
                command = lambda _: station_code_refresh())
        
        stationMenu.grid(row=1, column=2)
        
        stationEntry = Entry(channelFrame, width=40)
        stationEntry.state(["disabled"])
        stationEntry.grid(row=2, column=2)
        
        Label(master=channelFrame, text="Location (if set to 'Other', specify or leave blank)", justify="left").grid(row=3, column=1, sticky="w", padx=5, pady=2)
        
        locationMenu = OptionMenu(
                channelFrame, 
                locationVar, 
                location, 
                *list(location_codes + ["Other"]), 
                command = lambda _: location_code_refresh())
        
        locationMenu.grid(row=3, column=2)
        
        locationEntry = Entry(channelFrame, width=40)
        locationEntry.state(["disabled"])
        locationEntry.grid(row=4, column=2)
        
        Label(master=channelFrame, text="Channel (if set to 'Other', specify)", justify="left").grid(row=5, column=1, sticky="w", padx=5, pady=2)
        
        channelMenu = OptionMenu(
                channelFrame, 
                channelVar, 
                channel, 
                *list(channel_codes + ["Other"]),
                command = lambda _: channel_code_refresh())
        
        channelMenu.grid(row=5, column=2)
        
        channelEntry = Entry(channelFrame, width=40)
        channelEntry.state(["disabled"])
        channelEntry.grid(row=6, column=2)
        
        Label(master=channelFrame, text="Epoch (Deployment) Dates (if set to 'Other', specify)", justify="left").grid(row=7, column=1, sticky="w", padx=5, pady=2)
        
        deploymentDatesMenu = OptionMenu(
                channelFrame, 
                deploymentDatesVar, 
                deployment_times_for_entry[0], 
                *list(deployment_times_for_entry + ["Other"]),
                command = lambda _: deployment_times_refresh())
        
        deploymentDatesMenu.grid(row=7, column=2)
        
        deploymentDatesEntry = Entry(channelFrame, width=40)
        deploymentDatesEntry.state(["disabled"])
        deploymentDatesEntry.grid(row=8, column=2)
        
        #More channel parameters
        Label(master=channelFrame, text="Site Name", justify="left").grid(row=11, column=1, sticky="w", padx=5, pady=2)
        siteNameEntry = Entry(channelFrame, width=40)
        siteNameEntry.grid(row=11, column=2)
        
        Label(master=channelFrame, text="Latitude [decimal]", justify="left").grid(row=12, column=1, sticky="w", padx=5, pady=2)
        latitudeEntry = Entry(channelFrame, width=40)
        latitudeEntry.grid(row=12, column=2)
        
        Label(master=channelFrame, text="Longitude [decimal]", justify="left").grid(row=13, column=1, sticky="w", padx=5, pady=2)
        longitudeEntry = Entry(channelFrame, width=40)
        longitudeEntry.grid(row=13, column=2)
        
        Label(master=channelFrame, text="Elevation [masl]", justify="left").grid(row=14, column=1, sticky="w", padx=5, pady=2)
        elevationEntry = Entry(channelFrame, width=40)
        elevationEntry.grid(row=14, column=2)
        
        sitename = self.inv.select(
                station=stationVar.get(), 
                location=location)[0][0].site.name
        latitude = self.inv.select(
                station=stationVar.get(), 
                location=location)[0][0][0].latitude
        longitude = self.inv.select(
                station=stationVar.get(), 
                location=location)[0][0][0].longitude
        elevation = self.inv.select(
                station=stationVar.get(), 
                location=location)[0][0][0].elevation
        siteNameEntry.insert(END, str(sitename))
        latitudeEntry.insert(END, str(latitude))
        longitudeEntry.insert(END, str(longitude))
        elevationEntry.insert(END, str(elevation))
        
        Label(master=channelFrame, text="Sample Rate (>= 0) [Hz]", justify="left").grid(row=15, column=1, sticky="w", padx=5, pady=2)
        sampleRateEntry = Entry(channelFrame, width=40)
        sampleRateEntry.grid(row=15, column=2)
        
        Label(master=channelFrame, text="Azimuth (0 <= A < 360, from North) [degrees]", justify="left").grid(row=16, column=1, sticky="w", padx=5, pady=2)
        azimuthEntry = Entry(channelFrame, width=40)
        azimuthEntry.grid(row=16, column=2)
        
        Label(master=channelFrame, text="Depth (>= 0) [m]", justify="left").grid(row=17, column=1, sticky="w", padx=5, pady=2)
        depthEntry = Entry(channelFrame, width=40)
        depthEntry.grid(row=17, column=2)
        
        Label(master=channelFrame, text="Dip (-90 <= D <= 90, from Horizontal) [degrees]", justify="left").grid(row=18, column=1, sticky="w", padx=5, pady=2)
        dipEntry = Entry(channelFrame, width=40)
        dipEntry.grid(row=18, column=2)
        
        samplerate = self.inv.select(
                station=stationVar.get(), 
                location=location,
                channel=channel, 
                starttime=starttime, 
                endtime=endtime)[0][0][0].sample_rate
        azimuth = self.inv.select(
                station=stationVar.get(), 
                location=location,
                channel=channel, 
                starttime=starttime, 
                endtime=endtime)[0][0][0].azimuth
        depth = self.inv.select(
                station=stationVar.get(), 
                location=location,
                channel=channel, 
                starttime=starttime, 
                endtime=endtime)[0][0][0].depth
        dip = self.inv.select(
                station=stationVar.get(), 
                location=location,
                channel=channel, 
                starttime=starttime, 
                endtime=endtime)[0][0][0].dip
        sampleRateEntry.insert(END, str(samplerate))
        azimuthEntry.insert(END, str(azimuth))
        depthEntry.insert(END, str(depth))
        dipEntry.insert(END, str(dip))
        
        #response information
        responseFrame = LabelFrame(popup, text="Instrumental Response and Calibration", relief=RIDGE)
        responseFrame.grid(row=2, column=1, padx=3, pady=3, sticky="w")
        
        Label(responseFrame, text="Choose SENSOR response file: ", justify="left").grid(row=1, column=1, sticky="w", padx=5, pady=2)
        sensorResponseEntry = Entry(responseFrame, width=40)
        sensorResponseEntry.grid(row=1, column=2)
        sensorResponseButton = Button(responseFrame, text="Open", command=lambda: openResponseFile(sensorResponseEntry, _type="sensor"))
        sensorResponseButton.config(width=11)
        sensorResponseButton.grid(row=1, column=3, sticky="w")
        
        Label(responseFrame, text="Choose DATALOGGER response file: ", justify="left").grid(row=2, column=1, sticky="w", padx=5, pady=2)
        dataloggerResponseEntry = Entry(responseFrame, width=40)
        dataloggerResponseEntry.grid(row=2, column=2)
        dataloggerResponseButton = Button(responseFrame, text="Open", command=lambda: openResponseFile(dataloggerResponseEntry, _type="datalogger"))
        dataloggerResponseButton.config(width=11)
        dataloggerResponseButton.grid(row=2, column=3, sticky="w")
        
        Label(responseFrame, text="DATALOGGER response built in to SENSOR response file (e.g., for TitanEA): ", justify="left").grid(row=3, column=1, columnspan=2, sticky="w", padx=5, pady=2)
        dataloggerResponseVar = IntVar()
        dataloggerResponseVar.set(0)
        dataloggerResponseCheckbutton = Checkbutton(responseFrame, variable=dataloggerResponseVar, command=lambda: disableDataloggerDialogBox(dataloggerResponseEntry, dataloggerResponseButton))
        dataloggerResponseCheckbutton.grid(row=3, column=3, sticky="w")
        
        Label(responseFrame, text="Use only the nominal SENSOR response (i.e., already calibrated): ", justify="left").grid(row=4, column=1, columnspan=2, sticky="w", padx=5, pady=2)
        sensorCalibrationVar = IntVar()
        sensorCalibrationVar.set(0)
        sensorCalibrationCheckbutton = Checkbutton(responseFrame, variable=sensorCalibrationVar, command=lambda: noCalibration(_type="sensor"))
        sensorCalibrationCheckbutton.grid(row=4, column=3, sticky="w")
        
        Label(responseFrame, text="Use only the nominal DATALOGGER response (i.e., already calibrated): ", justify="left").grid(row=5, column=1, columnspan=2, sticky="w", padx=5, pady=2)
        dataloggerCalibrationVar = IntVar()
        dataloggerCalibrationVar.set(0)
        dataloggerCalibrationCheckbutton = Checkbutton(responseFrame, variable=dataloggerCalibrationVar, command=lambda: noCalibration(_type="datalogger"))
        dataloggerCalibrationCheckbutton.grid(row=5, column=3, sticky="w")
        
        #response information
        supportingMetadataFrame = LabelFrame(popup, text="Supporting Metadata", relief=RIDGE)
        supportingMetadataFrame.grid(row=3, column=1, padx=3, pady=3, sticky="w")
        
        Label(supportingMetadataFrame, text="Instrument type (e.g., Geophysical): ", justify = "left").grid(row=1, column=1, sticky="w", padx=5, pady=2)
        Label(supportingMetadataFrame, text="Sensor description (e.g., GeoSENSE BH-1 Corehole Seismometer): ", justify = "left").grid(row=2, column=1, sticky="w", padx=5, pady=2)
        Label(supportingMetadataFrame, text="Datalogger description (e.g., Guralp DM24-MK3 Datalogger): ", justify = "left").grid(row=3, column=1, sticky="w", padx=5, pady=2)
        Label(supportingMetadataFrame, text="Sensor Serial Number* (e.g., 116): ", justify = "left").grid(row=4, column=1, sticky="w", padx=5, pady=2)
        Label(supportingMetadataFrame, text="Datalogger Serial Number* (e.g., A1350): ", justify = "left").grid(row=5, column=1, sticky="w", padx=5, pady=2)
        Label(supportingMetadataFrame, text="Device ID*: ", justify = "left").grid(row=6, column=1, sticky="w", padx=5, pady=2)
        Label(supportingMetadataFrame, text="Website*: ", justify = "left").grid(row=7, column=1, sticky="w", padx=5, pady=2)

        instrumentTypeEntry = Entry(supportingMetadataFrame, width=30)
        instrumentTypeEntry.grid(row=1, column=2, columnspan=2)
        
        sensorDescriptionEntry = Entry(supportingMetadataFrame, width=30)
        sensorDescriptionEntry.grid(row=2, column=2, columnspan=2)
        
        dataloggerDescriptionEntry = Entry(supportingMetadataFrame, width=30)
        dataloggerDescriptionEntry.grid(row=3, column=2, columnspan=2)
        
        sensorSerialEntry = Entry(supportingMetadataFrame, width=30)
        sensorSerialEntry.grid(row=4, column=2, columnspan=2)
        
        dataloggerSerialEntry = Entry(supportingMetadataFrame, width=30)
        dataloggerSerialEntry.grid(row=5, column=2, columnspan=2)
        
        deviceIDEntry = Entry(supportingMetadataFrame, width=30)
        deviceIDEntry.grid(row=6, column=2, columnspan=2)
        
        websiteEntry = Entry(supportingMetadataFrame, width=30)
        websiteEntry.grid(row=7, column=2, columnspan=2)
        
        try:
            instrument_types = self.inv.select(
                    station=stationVar.get(), 
                    location=location, 
                    channel=channel, 
                    starttime=starttime, 
                    endtime=endtime)[0][0][0].types[0]
        except:
            instrument_types = ""
            
        try:
            sensor_description = self.inv.select(
                    station=stationVar.get(), 
                    location=location, 
                    channel=channel, 
                    starttime=starttime, 
                    endtime=endtime)[0][0][0].sensor.description.split("/")[0]
        except:
            sensor_description = ""
            
        try:
            datalogger_description = self.inv.select(
                    station=stationVar.get(), 
                    location=location, 
                    channel=channel, 
                    starttime=starttime, 
                    endtime=endtime)[0][0][0].sensor.description.split("/")[1]
        except:
            datalogger_description = ""
        
        try:
            sensor_serial_number = self.inv.select(
                    station=stationVar.get(), 
                    location=location, 
                    channel=channel, 
                    starttime=starttime, 
                    endtime=endtime)[0][0][0].sensor.serial_number.split("/")[0]
        except:
            sensor_serial_number = ""
        
        try:
            datalogger_serial_number = self.inv.select(
                    station=stationVar.get(), 
                    location=location, 
                    channel=channel, 
                    starttime=starttime, 
                    endtime=endtime)[0][0][0].sensor.serial_number.split("/")[1]
        except:
            datalogger_serial_number = ""
        
        try:
            external_references = self.inv.select(
                    station=stationVar.get(), 
                    location=location, 
                    channel=channel, 
                    starttime=starttime, 
                    endtime=endtime)[0][0][0].external_references
            if "Data Search" in external_references[0].description:
                website = external_references[0].uri
            else:
                website = external_references[1].uri
            if "Device" in external_references[0].description:
                device_id = external_references[0].uri.split("=")[1]
            else:
                device_id = external_references[1].uri.split("=")[1]
        except:
            device_id = ""
            website = ""
            
        instrumentTypeEntry.insert(END, str(instrument_types))
        sensorDescriptionEntry.insert(END, str(sensor_description))
        dataloggerDescriptionEntry.insert(END, str(datalogger_description))
        sensorSerialEntry.insert(END, str(sensor_serial_number))
        dataloggerSerialEntry.insert(END, str(datalogger_serial_number))
        deviceIDEntry.insert(END, str(device_id))
        websiteEntry.insert(END, str(website))

        Label(supportingMetadataFrame, text="* = optional entry", justify = "right").grid(row=100, column=3, sticky="e", padx=5, pady=2)
        
        #append
        appendMetaDataFrame = LabelFrame(popup, text="Submit Channel Metadata")
        appendMetaDataFrame.grid(row=4, column=1, padx=3, pady=3, sticky="w")
        
        appendMetaDataButton = Button(appendMetaDataFrame, text="Append to Inventory", command=lambda: appendToInventory())
        appendMetaDataButton.pack()
        
    def edit_inventory_entry(self):
        
        popup = Toplevel()
        popup.grab_set()
        popup.title("Edit Inventory Entry")
        
    def remove_inventory_entry(self):
        
        popup = Toplevel()
        popup.grab_set()
        popup.title("Remove Inventory Entry")

a = App("")
a.root.mainloop()