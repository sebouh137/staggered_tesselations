hcal_start_z=359.6

class Beampipe:
    def holeX(self, layer):
        z=hcal_start_z+layer*2.34
        return -1.329+(z-122.864)*(-11.252+1.329)/(522.664-122.864)
    def holeR(self, layer):
        z=hcal_start_z+layer*2.34
        clearance=3.85
        beampipe_thickness=0.2
        return clearance + beampipe_thickness + 4.403+(z-122.864)*(14.592-4.403)/(522.664-122.864)
    def beampipeX(self, layer):
        z=hcal_start_z+layer*2.34
        return -1.320+(z-122.864)*(-11.252+1.329)/(522.664-122.864)
    def beampipeR(self, layer):
        z=hcal_start_z+layer*2.34
        return 4.403+(z-122.864)*(14.592-4.403)/(522.664-122.864)
    def getZ(self, layer):
        return hcal_start_z+layer*2.34
