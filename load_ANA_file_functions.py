#  usefull conversions
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
if len(calnames)>2:
    print("Initializer::\tI've opened the files:\n",fnames[0].split('/')[-1],', ... ,',fnames[-1].split('/')[-1], '\t(total:',len(fnames),') and\n',
      calnames[0].split('/')[-1],', ... ,',calnames[-1].split('/')[-1], '\t(total:',len(calnames),')')
