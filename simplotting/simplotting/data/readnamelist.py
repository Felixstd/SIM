"""
Code used to read the namelist and assess the model parameters for 
the analysis. 

FSTD 15/05/2025
"""

class namelist:
    
    def __init__(self, configuration_exp = None, \
                       configuration_rheo = None, 
                       configuration_fig = None,
                       configuration_time = None):
        
        #--- Reading the Model Parameters ---#
        
        #------ Experiment ------#
        self.expno = [int(x) for x in configuration_exp['expno'].split(',')]
        self.savevar = int(configuration_exp['savevar'])
        self.plotfields = int(configuration_exp['plot_fields'])
        self.plotdilatation= int(configuration_exp['dilat'])
        self.read_all = int(configuration_exp['read_all'])
        self.plotsingle = int(configuration_exp['plot_single'])
        self.veltransect = int(configuration_exp['velocity_transect'])
        self.tracertransect = int(configuration_exp['tracers_transect'])
        
        #------ Rheology ------#
        self.muphi      = int(configuration_rheo['muphi'])
        self.mu0        = float(configuration_rheo['mu0'])
        self.muinf      = float(configuration_rheo['muinf'])
        self.mub        = [float(x) for x in configuration_rheo['mub'].split(',')]
        self.microangle = float(configuration_rheo['micro_angle'])
        self.dx         = float(configuration_rheo['dx'])*1e3
        self.Nx         = float(configuration_rheo['Nx'])
        self.Ny         = float(configuration_rheo['Ny'])
        
        #------ Figures ------#
        self.outputdir = str(configuration_fig['outputdir'])
        self.figdir    = str(configuration_fig['figdir'])
        
        #------ Time -------#
        self.dt = float(configuration_time['Dtminute'])
        self.startk = int(configuration_time['start_k'])
        
        