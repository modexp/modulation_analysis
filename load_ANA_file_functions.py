#  usefull conversions
import ROOT                        # To open root files
import numpy as np                 # Fast
import matplotlib.pyplot as plt    # Plotting
import matplotlib.cm as cm         # colorbars
import matplotlib                  # Plotting
import datetime                    # For dates on axis
from os import listdir             # To see what is in a folder
from scipy.optimize import curve_fit # Fitting
from scipy import interpolate
from scipy.optimize import optimize # Fitting
from scipy import signal
import emcee
import corner
import scipy.optimize as op
from matplotlib.ticker import MaxNLocator
import io
from iminuit import Minuit, describe, Struct

y2s      = 365.25 * 24 * 3600 # years to seconds
s2y      = 1. / y2s         # seconds to years
toffset  = 70 * y2s       # the timestamp from labview is 70 years off
droppath = '/home/tmons/Documents/masterthesis/plots/' # this is where I want to save my plots

matplotlib.rc('font', size=24)                   # Use big fonts...
plt.rcParams['figure.figsize'] = (12.0, 10.0)    # ... and big plots

# Dictonary peak format str(channel) + str(peakid) = E (keV)
pkE = {
    '00':1460,
    '10':1460,
    '20':511,   '21':1157.020,  '22':511+1157.020,
    '30':511,   '31':1157.020,  '32':511+1157.020,
    '40':1173.2,'41':1332.5,    '42':1173.2+1332.5,
    '50':1173.2,'51':1332.5,    '52':1173.2+1332.5,
    '60':661.7,
    '70':661.7}
# A list of the energy of the peak in the spectrum sourceE[channel][peakid] = E (keV)
sourceE = [[1460],
           [1460],
           [511,   1157.020, 511 + 1157.020],
           [511,   1157.020, 511 + 1157.020],
           [1173.2,1332.5,   1173.2 + 1332.5],
           [1173.2,1332.5,   1173.2 + 1332.5],
           [661.7],
           [661.7]]
# A list of the halflife of the source format halflife[channel] = tau1/2
halflife = [1e9,    1e9,
            59.1,   59.1,
            5.2711, 5.2711,
            30.08,  30.08]
dhalflife = [1e9,   1e9,
            0.3,    0.3,
            0.0004, 0.0004,
            0.09,   0.09]

sourceName = ['Background',
           'Background',
            'Ti-44',
            'Ti-44',
            'Co-60',
            'Co-60',
            'Cs-137',
            'Cs-137']

sourceLatex =['Background',
           'Background',
            '$^{44}$Ti',
            '$^{44}$Ti',
            '$^{60}$Co',
            '$^{60}$Co',
            '$^{137}$Cs',
            '$^{137}$Cs']


col_lis    = ['green', 'red', 'blue', 'cyan', 'brown', 'pink', 'magenta','black']

all_files = sorted(listdir(datafolder)) # All files in this directory

physics_exclude = [
    # Rates in all channel too small. Might be a problem in delta_t or some downtime
    'ANA_mx_n_20160817_1855.root',
    # Problem with channel 7. Calibration constants changing too rapidly and we get strange rates henceforth
    'ANA_mx_n_20160530_0926.root', 'ANA_mx_n_20160601_0811.root',
    # Problem with Ti-channels, the high end of the energy spectrum is off
    'ANA_mx_n_20151130_0922.root', 'ANA_mx_n_20151202_0906.root', 'ANA_mx_n_20151204_1911.root',
    'ANA_mx_n_20151207_0819.root',
    # Double rates???
    'ANA_mx_n_20161028_0726.root',
    # The HV measurements
    'ANA_mx_n_20170220_0945.root', 'ANA_mx_n_20170220_1130.root', 'ANA_mx_n_20170220_1246.root',
    'ANA_mx_n_20170220_1444.root', 'ANA_mx_n_20170220_1900.root', 'ANA_mx_n_20170220_2027.root',
    'ANA_mx_n_20170221_0714.root', 'ANA_mx_n_20170221_0848.root', 'ANA_mx_n_20170221_1101.root',
    'ANA_mx_n_20170221_1242.root', 'ANA_mx_n_20170221_1426.root', 'ANA_mx_n_20170221_1607.root',
    'CAL_mx_n_20170220_0945.root', 'CAL_mx_n_20170220_1130.root', 'CAL_mx_n_20170220_1246.root',
    'CAL_mx_n_20170220_1444.root', 'CAL_mx_n_20170220_1900.root', 'CAL_mx_n_20170220_2027.root',
    'CAL_mx_n_20170221_0714.root', 'CAL_mx_n_20170221_0848.root', 'CAL_mx_n_20170221_1101.root',
    'CAL_mx_n_20170221_1242.root', 'CAL_mx_n_20170221_1426.root', 'CAL_mx_n_20170221_1607.root',
    'ANA_mx_n_20170307_1346.root', 'ANA_mx_n_20170307_1550.root', 'ANA_mx_n_20170307_1551.root',
    'ANA_mx_n_20170307_1753.root', 'ANA_mx_n_20170307_2037.root', 'ANA_mx_n_20170308_0709.root',
    'ANA_mx_n_20170308_0930.root', 'ANA_mx_n_20170308_1152.root', 'ANA_mx_n_20170308_1401.root',
    'ANA_mx_n_20170308_1600.root', 'ANA_mx_n_20170308_1759.root',
#     'ANA_mx_n_20170308_1759.root',
#     'CAL_mx_n_20170306_0823.root',
    'CAL_mx_n_20170307_1346.root', 'CAL_mx_n_20170307_1550.root', 'CAL_mx_n_20170307_1551.root',
    'CAL_mx_n_20170307_1753.root', 'CAL_mx_n_20170307_2037.root', 'CAL_mx_n_20170308_0709.root',
    'CAL_mx_n_20170308_0930.root', 'CAL_mx_n_20170308_1152.root', 'CAL_mx_n_20170308_1401.root',
    'CAL_mx_n_20170308_1600.root', 'CAL_mx_n_20170308_1759.root',
    # disk space full, problems with writing files to disk
    'ANA_mx_n_20170413_0720.root', 'CAL_mx_n_20170413_0720.root',
#     'ANA_mx_n_20170411_1033.root', 'CAL_mx_n_20170411_1033.root'
    'ANA_mx_n_20170421_0747.root', 'ANA_mx_n_20170508_0709.root', 'CAL_mx_n_20170508_0709.root',
# 'ANA_mx_n_20170508_0709.root',
    'CAL_mx_n_20170509_0743.root', 'ANA_mx_n_20170510_1503.root', 'ANA_mx_n_20150701_1303.root',
    'ANA_mx_n_20150702_1150.root', 'ANA_mx_n_20150702_1349.root', 'ANA_mx_n_20150707_1248.root',
    'ANA_mx_n_20150708_0755.root', 'ANA_mx_n_20150708_0947.root', 'ANA_mx_n_20150709_0727.root',
    'ANA_mx_n_20160617_0958.root', # All rates get to small in this file
    'ANA_mx_n_20160616_0951.root',  # All rates get to small in this file
    'ANA_mx_n_20161122_0758.root', 'ANA_mx_n_20161121_1433.root', 'ANA_mx_n_20161123_1218.root']

# The files that did not work, trying to open thise files will result in a problem later on. There is
badlist = [#'ANA_mx_n_20151016_1419.root',
           'ANA_mx_n_20160412_0739.root', 'ANA_mx_n_20160418_0950.root', 'ANA_mx_n_20160527_1421.root',
           'CAL_mx_n_20151016_1419.root', 'ANA_mx_n_20151016_1419.root', 'CAL_mx_n_20151022_1134.root',
           'CAL_mx_n_20151022_1143.root', 'CAL_mx_n_20151022_1200.root', 'CAL_mx_n_20151022_1207.root',
           'CAL_mx_n_20160412_0739.root', 'CAL_mx_n_20160418_0950.root', 'CAL_mx_n_20160527_1421.root',
           'CAL_mx_n_20161215_1047.root', 'ANA_mx_n_20161215_1047.root', 'ANA_mx_n_20170105_0836.root',
           'CAL_mx_n_20170105_0836.root',
          # We had a power failure therefore exlude (there is also no data so it is not interesting:
          'CAL_mx_n_20170117_1607.root', 'CAL_mx_n_20170117_1620.root', 'CAL_mx_n_20170117_1648.root',
          'CAL_mx_n_20170117_1653.root', 'CAL_mx_n_20170117_1942.root', 'CAL_mx_n_20170118_1634.root',
          # An I/O error caused LabView to crash
          'CAL_mx_n_20170127_1007.root']

new_physics_exclude= ['ANA_mx_n_20170629_0710.root', 'ANA_mx_n_20170629_1140.root', 'ANA_mx_n_20171017_1144.root',
                      'ANA_mx_n_20171027_1201.root', 'ANA_mx_n_20180115_1258.root', 'ANA_mx_n_20180116_1142.root',
                      'ANA_mx_n_20180117_0835.root', 'ANA_mx_n_20180119_1215.root', 'ANA_mx_n_20180122_1503.root',
                      'ANA_mx_n_20180125_1018.root', 'ANA_mx_n_20180126_1218.root', 'ANA_mx_n_20180219_0914.root',
                      'ANA_mx_n_20180220_1148.root', 'ANA_mx_n_20180221_1042.root', 'ANA_mx_n_20180222_1031.root',
                      'ANA_mx_n_20180913_0747.root', 'ANA_mx_n_20180913_1059.root', 'ANA_mx_n_20180914_0704.root',
                      'CAL_mx_n_20180914_0704.root']

new_badruns_exclude= []

badlist += physics_exclude
badlist += new_physics_exclude

def files_to_open(set_date, i_want='ANA'):
    '''Opens files in this directory. I have CAL files and ANA files. I start from the first date as given
    by the first argument'''
    fnames=[]

    # Loop through all the files
    for file in all_files:
        if file=='.directory':
            continue
        try:
            test=int(file[9:15])
        except:
            continue
        if int(file[9:15])>=fromdate and int(file[9:15])<=todate:
            if file[0:3]==i_want and file not in badlist:
                if (int(file[9:17])>=20170511 and int(file[9:17])<=20170630) or (int(file[9:17])>=20160501 and int(file[9:17])<=20160915):
                    continue
                fnames.append(datafolder + file)
    return fnames

fnames = files_to_open(fromdate)
calnames = files_to_open(fromdate,i_want='CAL')
#if len(calnames)>2:
#    print("Initializer::\tI've opened the files:\n",fnames[0].split('/')[-1],', ... ,',fnames[-1].split('/')[-1], '\t(total:',len(fnames),') and\n',
#      calnames[0].split('/')[-1],', ... ,',calnames[-1].split('/')[-1], '\t(total:',len(calnames),')')


in_ana = ['t0', 'time', 'channel', 'peak', 'rate', 'drate', 'e', 'res', 'temp', 'pres', 'bx', 'by',
          'bz', 'btot', 'humid','frac','hv0', 'hv1','hv2','hv3','hv4', 'hv5','hv6','hv7','chi2ndf',
          'bgrate', 'bgdrate', 'tot_error']
in_cal = ['id', 'cal_tmin', 'cal_tmax', 'c0', 'c1', 'c2', 'chi2' ]

time_start = []
time_end   = []

def file_opener(fnames = files_to_open('201509'), calnames =  files_to_open('201509','CAL')):
    badfiles = []
    for k in in_ana: # For all properties from in_ana make a list of this property
        globals()[k] = []

    # For each file in fnames open it and append the property to the right list
    for i in range(len(fnames)):
        f = ROOT.TFile.Open(fnames[i]) # Open the ANA file
        try: # Sometimes the anafiles are corrupted, there will be no ana;1.
            tree = f.Get("ana;1")
        except ReferenceError: # If this happens, a reference error occurs, add this name to a list of bad files
            badfiles.append(fnames[i])
            continue
        # A file can also be corrupted such that it is a TObject instead of a tree with leaves
        if 'TObject' in str(type(tree)):
            badfiles.append(fnames[i])
            continue
        # For each event in the tree, get its properties and append to the right list
        for j , event in enumerate( tree) :
            for k in in_ana:
                if k == 'time': # Time is saved wrt t0, therefore need to be added to get the absolute time
                    if j == 0: dt = (getattr(event, k))
                    timestamp = getattr(event, k) + getattr(event, 't0')
                    globals()[k].append(timestamp)
                    time_start.append(timestamp - dt)
                    time_end.append(  timestamp + dt)
                elif k == 'frac': # for frac some events do not have this property, add a value of None here.
                    try:
                        globals()[k].append( getattr(event, k) )
                    except AttributeError:
                        globals()[k].append( None )
                elif k == 'bgrate' or k == 'bgdrate': # for chi2ndf some events do not have this property, add a value of None here.
                    try:
                        globals()[k].append( getattr(event, k) )
                    except AttributeError:
                        globals()[k].append( -1 )
                elif k == 'chi2ndf': # for chi2ndf some events do not have this property, add a value of None here.
                    try:
                        globals()[k].append( getattr(event, k) )
                    except AttributeError:
                        globals()[k].append( -1 )
                elif k == 'tot_error': # for tot_error some events do not have this property, add a value of None here.
                    try:
                        globals()[k].append( getattr(event, k) )
                    except AttributeError:
                        globals()[k].append( -1 )
                elif 'hv' in k: # for hv0-hv7 some events do not have this property, add a value of None here.
                        try:
                            globals()[k].append( getattr(event, k) )
                        except AttributeError:
                            globals()[k].append( None )
                elif 'drate' in k:
                    if getattr(event, k) == 'nan' or getattr(event, k) == 'NaN':
                        print('Initializer::file_opener::A "None" occured')
                        globals()[k].append( -10 )
                    else:
                        globals()[k].append( getattr(event, k) )
                else: # Most poperties get added here to the right list
                    globals()[k].append( getattr(event, k) )


    for k in in_cal: # For all properties from in_cal make a list of this property
        globals()[k] = []

    # For each file in calnames open it and append the property to the right list
    for i in range(len(calnames)):
        f = ROOT.TFile.Open(calnames[i]) # open the CAL file
        try: # Sometimes the calibration files are corrupted, there will be no cal;1.
            tree = f.Get("cal;1")
        except ReferenceError: #If this happens, a reference error occurs, add this name to a list of bad files
                badfiles.append(calnames[i])
                continue
        # For each event in the tree, get its properties and append to the right list
        for j , event in enumerate(tree) :
            for k in in_cal:
                # Some of the properties are not a simple leaf in the (cal) root file but are root vectors
                # these are especially tricky to open. These can be opened by treating the rootvector as a
                # variable and than converting it to a list before adding it to the list of this property
                if k == 'c0' or k == 'c1' or k == 'c2' or k =='chi2':
                    vec1 = getattr(event, k)
                    list1 = [value for value in vec1] # This
                    globals()[k].append(list1)
                else:
                    # Sometimes strange things happend with the calibration time therefore I check if the
                    # timestamp does not have strange values.
                    if k == 'cal_tmin':
                        if getattr(event, k) < 1:
                            print('cal time is negative')
                            badfiles.append(calnames[i])
                        elif getattr(event, k) > 4e9:
                            print('cal time to big in file %s' %calnames[i])
                            if calnames[i] not in badfiles:
                                badfiles.append(calnames[i])
                        else:
                            globals()[k].append( getattr(event, k) )
    # As written below, we now have some properties from the calibration that are in a vector. To convert these
    # to a list we first create some lists:
    for string in ['c0','c1', 'c2', 'chi2']:
        for j in range(8):
            name = string + '_chan_' + str(j)
            globals()[name] = []
            for i in range(len(c0)):
                ccc = globals()[string]
                globals()[name].append(ccc[i][j])

    # And then convert these to np arrays.
    for cx in ['c0','c1', 'c2', 'chi2']:
        for ch in range(8):
            name = cx + '_chan_' + str(ch)
            in_cal.append(name)

    for k in in_cal:
        if k == 'id':
            id_array = np.array(globals()[k])
        else:
            globals()[k] = np.array(globals()[k])
    in_ana.append('time_start')
    in_ana.append('time_end')
    for k in in_ana:
        globals()[k] = np.array(globals()[k])

    if not badfiles == []: print ('Initializer::The files that have a problem are:\t', badfiles)
    print ('Initializer::\tLoading done\n\tWe can access the properties:\n\t', *in_ana + in_cal)
    print ('Initializer::\tEach of these properties is an np.array and can be used accourdingly')
#file_opener(fnames, calnames)
