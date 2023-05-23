#all parameters in the project
import numpy as np
from shapely import *


sidelengths=[1.889]*9+[3.589]*9+[4.486]*36

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
class BeamPipe():
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
beampipe0=BeamPipe()

    
#determine the layer boundaries as a polygon
def layer_boundaries(layer=0, side="L", beampipe=beampipe0, height=det_height, width=det_width,gap_between_sides=0.4,
                    gap_for_backplanes=0.4, backplanes_in_holes=False):
    gap=gap_between_sides+gap_for_backplanes
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