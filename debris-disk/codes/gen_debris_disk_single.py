"""
Generate debris disk
@author: Eve J. Lee
Feb. 16th 2016
"""

import numpy as np
import debris_disk_single
import consts
import sys

#response = input("File name: ")
response = sys.argv[1]

dd = debris_disk_single.DebrisDisk("../disk_input/%s.txt"%response)
dd.AddSinglePlanet()
#dd.ComputeParentSingle(manual=True, amin=62.6, amax=75.9, I0=0, e0=0.)
dd.ComputeParentSingle(Nparticles=1)
dd.ComputeParentOrbital()
#dd.ComputeDustGrains()
dd.ComputeDustGrains_BetaOptimized()
dd.OutputParentOrbit("../parentorbit/%s"%response)

dd.ComputeBackgroundParentSingle()
dd.ComputeBackgroundParentOrbital()
dd.ComputeBackgroundDustGrains()

#dd.Plot_pqhk()
#dd.Plot3DOrbit()


dd.OutputDustAndBackOrbit("../dustorbit/%s"%response)













# """
# Generate debris disk
#
# @author: Eve J. Lee
# Feb. 16th 2016
# """
#
# import numpy as np
# import debris_disk_single
# import consts
# import sys
#
# #response = input("File name: ")
# response = sys.argv[1]
#
# dd = debris_disk_single.DebrisDisk("../disk_input/%s.txt"%response)
# dd.AddSinglePlanet()
# #dd.ComputeParentSingle(manual=True, amin=62.6, amax=75.9, I0=0, e0=0.)
# dd.ComputeParentSingle(Nparticles=1)
# dd.ComputeParentOrbital()
# dd.ComputeDustGrains()
#
# #dd.Plot_pqhk()
# #dd.Plot3DOrbit()
#
# dd.OutputParentOrbit("../parentorbit/%s"%response)
# dd.OutputDustOrbit("../dustorbit/%s"%response)
#