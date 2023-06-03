#all parameters in the project
import numpy as np
from shapely import *

#not including the last absorber
n_layers=64


nsmall=7;nmed=9;nlarge=n_layers-nsmall-nmed
sidelengths=[1.889]*nsmall+[3.589]*nmed+[4.486]*nlarge

#= 19 small cells tall + thickness of a cell wall
# or 10 medium cells tall + thickness of a cell wall
# or 8 large cells tall + thickness of a cell wall.  
det_height=62.24 

#12 inches
det_width=60.96

#gap between left and right sides of absorbers (this space includes the additional space for backplanes
det_gap=1.2

#first layer of the first absorber
hcal_start_z=359.6



class Component():
    def __init__(self,z_offset, thickness, name):
        self.z_offset=z_offset
        self.thickness=thickness
        self.name=name
z=0
components={}
#insulator may be slightly thinner in reality (66 microns rather than 100 for one possible product); in which case the air gaps
#would be slightly larger
for name, thickness in [("absorber", 1.52), ("airgap1", 0.0365),("cover", 0.04),                     
                        ("foil1", 0.015),("scint", 0.3), ("foil2", 0.015), ("pcb", 0.08),
                        ("insulator", 0.007),
                        ("airgap2", 0.0365)]:
    components[name]=Component(z, thickness, name)
    z+=thickness

coord_thickness=z# reduced to 2.11 due to thinner PCB and cover # NIM paper had 2.34
#del z

class BeamPipe():
    def holeX(self, layer):
        z=hcal_start_z+layer*coord_thickness
        return -1.329+(z-122.864)*(-11.252+1.329)/(522.664-122.864)
    def holeR(self, layer):
        z=hcal_start_z+layer*coord_thickness
        clearance=3.85
        beampipe_thickness=0.2
        return clearance + beampipe_thickness + 4.403+(z-122.864)*(14.592-4.403)/(522.664-122.864)
    def beampipeX(self, layer):
        z=hcal_start_z+layer*coord_thickness
        return -1.320+(z-122.864)*(-11.252+1.329)/(522.664-122.864)
    def beampipeR(self, layer):
        z=hcal_start_z+layer*coord_thickness
        return 4.403+(z-122.864)*(14.592-4.403)/(522.664-122.864)
    def getZ(self, layer):
        return hcal_start_z+layer*coord_thickness
beampipe0=BeamPipe()

    
#determine the layer boundaries as a polygon
def layer_boundaries(layer=0, side="L", beampipe=beampipe0, height=det_height, width=det_width,gap_between_sides=0.4,
                    gap_for_backplanes=0.4, backplanes_in_holes=False):
    gap=gap_between_sides+2*gap_for_backplanes
    offsetX=-10
    holeX=beampipe.holeX(layer)
    holeR=beampipe.holeR(layer)
    if side=="R":
        offsetX-width/2
        
        phi= np.linspace(-np.pi/2, np.pi/2, 25)
        
        x,y = [offsetX-width/2, offsetX-width/2, -gap/2, -gap/2] + list(holeX-holeR*np.cos(phi)) + [-gap/2,-gap/2, -39.8], \
                 [-height/2,height/2, height/2, holeR]+ list(-holeR*np.sin(phi))+ [-holeR,-height/2, -height/2]
        poly = Polygon(zip(x,y))
        if backplanes_in_holes:
            #remove some space for the backplanes in the sides of the detector:
            x=[0, holeX, holeX, 0, 0]
            y=holeR+gap_for_backplanes
            y=[y,y, -y, -y, y]
            return poly.difference(Polygon(zip(x,y)))
        else :
            return poly
    if side=="L":
        
        phi0 = np.arccos((holeX-gap/2)/holeR)
        #print((np.pi-phi0)/np.pi*2)
        phi = np.linspace(phi0,2*np.pi-phi0, 13)
        
        x=[offsetX+width/2, offsetX+width/2, gap/2] + list(holeX-np.cos(phi)*holeR) + [gap/2, 19.8]
        y=[-height/2,height/2, height/2]+list(holeR*np.sin(phi)) + [-height/2, -height/2]
        return Polygon(zip(x,y))
    
#thickness of the cell walls in the 3d-printed frame
wall_thickness = 0.08
#gap between the cell cell of the 3d printed frame and the scintillator tile.  
wall_scint_gap=0.01

nW=34
absorber_material=['W']*nW+['Fe']*(n_layers-nW+1)