import debris_disk_single
import sys

response = sys.argv[1]

dd = debris_disk_single.DebrisDisk("../disk_input/%s.txt"%response)

dd.AddSinglePlanet()
dd.ComputeParentSingle(Nparticles=1, coll_in_middle=True)
dd.ComputeParentOrbital()
dd.ComputeDustGrains_Optimized()
dd.OutputParentOrbit("../parentorbit/%s"%response)

dd.ComputeBackgroundParentSingle_Optimized()
dd.ComputeBackgroundParentOrbital_Optimized()
dd.ComputeBackgroundDustGrains_Optimized()

dd.ComputeForkDust_Optimized3(response)

dd.OutputDustAndBackOrbit("../dustorbit/%s"%response)

