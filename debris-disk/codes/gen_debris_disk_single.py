import debris_disk_single
import sys

response = sys.argv[1]

dd = debris_disk_single.DebrisDisk("../disk_input/%s.txt"%response)

dd.AddSinglePlanet()
dd.ComputeParentSingle(Nparticles=1, coll_in_middle=True)
dd.ComputeParentOrbital()
dd.ComputeDustGrains()
dd.OutputParentOrbit("../parentorbit/%s"%response)

dd.ComputeBackgroundParentSingle()
dd.ComputeBackgroundParentOrbital()
dd.ComputeBackgroundDustGrains()

dd.ComputeForkDust()

dd.OutputDustAndBackOrbit("../dustorbit/%s"%response)

