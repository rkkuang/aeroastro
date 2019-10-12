from utils import gen_site

#http://vlbiimaging.csail.mit.edu/
#This form allows you to generate your own simulated VLBI data from specified array and target parameters and your own source image! Use this to debug your algorithm or to generate realistic data similar to your problem. Once you submit the form we will generate an OIFITS file with the simulated visiblities and bispectrum data along with a number of plots for visualization. Simulated data is generated using the MAPS software package.

#Step 1: Select Image of the Emission
#Select or upload an image that you would like to observe and specify the total flux density of the emission.
Total_flux_density = 1 #Janskys
Rotation = 0 #Degrees

Source_image = []#SgrA* Model
Spin = 0% # 0-100 percentage
Inclination = 0# 0-89 degree


#Step 2: Select Direction and FOV
#Identify the direction to the target source. Right ascension should be in the form HH:MM:SS.SS for hours, minutes, and seconds and declination should be in the form DD:MM:SS.SS for degrees, arcminutes, and arcseconds. Field of view is specified in arcseconds. Warning: You must choose coordinates such that your region will be observable from your observatory site (the first telescope you specify below) at the start time that you specify, otherwise the resulting output will be incorrect.

#Field Of View Center: Right Ascension (HH:MM:SS.SS) 17:45:40.041 Declination(DD:MM:SS.SS) -29:00:28.118
FideldofViewCenter_RightAscesion = "17:45:40.041"
FideldofViewCenter_Declination = "-29:00:28.118"

#Field Of View Size:         Right Ascension (arcseconds) 0.000204595   Declination (arcseconds) 0.000204595
FieldofViewSize_RightAscesion = 0.000204595
FieldofViewSize_Delination = 0.000204595

'''
Step 3: Specify Telescope Array
Add the telescope locations and intrinsic parameters that you would like to use to simulate data

Initilization(typo): Select a pre-loaded telescope
Name: Unique name for each telescope station (up to 12 characters)
East Longitude/Latitude: East longitude and latitude of the array center. For locations less than 180 degrees west of Greenwich a minus sign should precede the longitude entry.
X/Y/Z Position: Absolute X, Y, Z coordinates of each station (in meters) relative to the center of the Earth
Lower/Upper Elevation: Lower and upper elevation limits of the of the antenna in degrees
SEFD: System equivalent flux denisty of the antenna
Diameter: Antenna diameter in meters
'''
# site = [Name, East Longitude,	Latitude, X-Position, Y-Position, Z-Position, Lower Elevation, Upper Elevation, SEFD, Diameter]
# X = R_E*cos(lati)*cos(longi)
# Y = R_E*cos(lati)*sin(longi)
# Z = R_E*sin(lati)
#LMT = {"name":"LMT","East_Longitude":"-97:18:53","Latitude":"18:59:06","X_position":-768713.9637,"Y_position":-5988541.7982,"Z_position":2063275.9472,"Lower_Elevation":15,"Upper_Elevation":85,"SEFD":560,"Diameter":50}
# gen_site(name, East_Longitude, Latitude , X_position, Y_position, Z_position , Lower_Elevation, Upper_Elevation, SEFD, Diameter)
LMT = gen_site("LMT", "-97:18:53", "18:59:06", -768713.9637, -5988541.7982, 2063275.9472, 15, 85, 560, 50)
PV = gen_site("PV", "-3:23:33.8", "37:03:58.2", 5088967.9, -301681.6, 3825015.8, 15, 85, 2900, 30)
ALMA50 = gen_site("ALMA50", "-67:45:11.4", "-23:01:09.4", 2225037.1851, -5441199.162, -2479303.4629, 15, 85, 110, 84.7)
#SMTO = gen_site("SMTO", "-109:52:19", "32:42:06", -1828796.2, -5054406.8, 3427865.2, 15, 85, 11900, 10)
SMT = gen_site("SMT", "-109:52:19", "32:42:06", -1828796.2, -5054406.8, 3427865.2, 15, 85, 11000, 10)
Hawaii8 = gen_site("Hawaii8", "-155:28:40.7", "19:49:27.4", -5464523.4, -2493147.08, 2150611.75, 15, 85, 4900, 20.8)
#PdBI = gen_site("PdBI", "05:54:28.5", "44:38:02.0", 4523998.4, 468045.24, 4460309.76, 15, 85, 1600, 36,7)
PDB = gen_site("PDB", "05:54:28.5", "44:38:02.0", 4523998.4, 468045.24, 4460309.76, 15, 85, 5200, 36,7)
SPT = gen_site("SPT", "-000:00:00.0", "-90:00:00", 0, 0, -6359587.3, 15, 85, 7300, 12)
GLT = gen_site("GLT", "72:35:46.4", "38:25:19.1", 1500692, -1191735, 6066409, 15, 85, 4744, 12)
CARMA8 = gen_site("CARMA8", "-118:08:30.3", "37:16:49.6", -2397431.3, -4482018.9, 3843524.5, 15, 85, 3500, 26.9)
SMA = gen_site("SMA", "-155:28:40.7", "19:49:27.4", -5464523.4, -2493147.08, 2150611.75, 15, 85, 4000, 20.8)

Array = [LMT,PV,ALMA50,SMT,Hawaii8,PDB,SPT,GLT,CARMA8,SMA]

'''
Step 4: Specify Date and Time Data is Collected
Specify the time of when you would like measurments to be taken, and the time interval between measurements.

Start Time: Specify the time of your first observation in Universal Time (UT). The required format is "YYYY:ddd:hh:mm:ss" where YYYY is the year, ddd is the day number (e.g., December 31 is day 365); hh is the UT hour, mm is the UT minute, and ss is the UT second.
Scan Duration: The length of a continuous scan in seconds
Interval Length: The time in seconds between successive scans
Number of Samples: The number of successive scans of this type
'''
#Start Time (UT)    Scan Duration (seconds) Interval Length (seconds)   Number of Samples
data_colle1 = ("2016:95:00:00:00", 10, 600, 100)
data_colle2 = ("2016:97:00:00:00", 10, 600, 100)

DATA_COLLE = [data_colle1, data_colle2]

'''
Step 5: Specify Collection Parameters
Specify the center frequency and width of the observing channel in MHz.
Center Frequency (MHz): 227297   Bandwidth (MHz): 4096

Specify your integration time in seconds (sometimes referred to as “dump time” or “record length”). This is not the total duration of your observation, but rather the sampling and recording interval of the data.

Integration Time (seconds): 10
'''
CenterFrequency = 227297 # MHz
Bandwidth = 4096 # MHz
IntegrationTime = 10 # seconds


'''
Step 6: Add Noise and Generate Data
 Simulate Without ANY Noise

 Simulate Without Atmospheric Phase Errors

 Simulate Without 5% Gain Error
'''
SimuWithoutAnyNoise = False 
SimuWithoutAtmosphericPhaseErrors = False
GainError = 0.05
SimuWithoutGainError = False