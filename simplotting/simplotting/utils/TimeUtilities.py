"""
Time Utility file made by Felix St-Denis to process 
model outputs.

The code was inspired by the time utility made by Mathieu Plante. 

"""

from datetime import datetime, date, time, timedelta

class TimeUtility:
    
    def __init__(self, configuration = None):
        
        if configuration != None:
            
            self.StartYear    = int(configuration['start_year'])
            self.StartMonth   = int(configuration['start_month'])
            self.StartDay     = int(configuration['start_day'])
            self.StartHour    = int(configuration['start_hour'])
            self.StartMinute  = int(configuration['start_minute'])
            self.EndYear      = int(configuration['end_year'])
            self.EndMonth     = int(configuration['end_month'])
            self.EndDay       = int(configuration['end_day'])
            self.EndHour      = int(configuration['end_hour'])
            self.EndMinute    = int(configuration['end_minute'])
            self.dt           = int(configuration['Dtminute'])
            
            self.StartDate = datetime(self.StartYear,self.StartMonth,self.StartDay, hour = self.StartHour, minute = self.StartMinute )
            self.EndDate   = datetime(self.EndYear,self.StartMonth,self.EndDay, hour = self.EndHour, minute = self.EndMinute )
            
            self.ndays     = int((self.EndDate - self.StartDate).days +1)
            self.nsteps    = int(int((self.EndDate - self.StartDate).days)*24.0*60/self.dt + int((self.EndDate - self.StartDate).seconds/(60*self.dt) +1))
            
    
    def dates_analysis(self, nmax = None):
            
        for n in range(0, nmax):
            yield self.StartDate + timedelta(minutes = n*self.dt)