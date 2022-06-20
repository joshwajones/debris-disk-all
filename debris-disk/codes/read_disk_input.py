"""
Read in debris disk model input

@author: Eve J. Lee
Feb. 17th 2016
"""

def ReadInput(fname):
    inputdata = {}
    with open(fname, 'r') as fileobject:
        for line in fileobject:
            label, data = line.split('=')
            inputdata[label.strip()]=float(data.split('\\')[0])
    return inputdata

def WriteInput(data,fname):
    f = open(fname, 'w')
    for key in data:
        f.write("%s = %f\n"%(key,data[key]))
    f.close()

    
