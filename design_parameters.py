#all parameters in the project
import numpy as np
from shapely import *

#not including the last absorber
n_layers=60


#nsmall=7;nmed=7;nlarge=n_layers-nsmall-nmed
#sidelengths=[1.889]*nsmall+[3.589]*nmed+[4.486]*nlarge

#test
nsmall=7;nlarge=0;nmed=n_layers-nlarge-nsmall; 
sidelengths=[1.889]*nsmall+[3.112]*nmed



#from the blueprints
det_height=59.29 

det_right_width=39.44

det_left_width=19.44

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
def layer_boundaries(layer=0, side="L", beampipe=beampipe0, backplanes_in_holes=False):
    
    gap=0.38#gap between left and right absorbers
    
    holeR=13.7+(17.04-13.7)*(layer-1)/64
    if side=="R":
        #relative to the inner edge of the absorber
        holeX=7.07+(10.33-7.07)*(layer-1)/64
        
        phi= np.linspace(-np.pi/2, np.pi/2, 25)
        
        x = np.concatenate([[-det_right_width, -det_right_width, 0, 0, -0.5, -0.5,0, 0],-holeX-holeR*np.cos(phi),[0,0,-0.5,-0.5,0,0, -det_right_width]])-gap/2
        y= np.concatenate([[-det_height/2,det_height/2, det_height/2,28.14, 28.14,18.14,18.14, holeR],-holeR*np.sin(phi),[-holeR,-18.14, -18.14,-28.14, -28.14,-det_height/2, -det_height/2]])
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
        #relative to the inner edge of the absorber
        holeX=7.44+(10.71-7.44)*(layer-1)/64
        
        phi0 = np.arccos(-holeX/holeR)
        #print((np.pi-phi0)/np.pi*2)
        phi = np.linspace(phi0,2*np.pi-phi0, 13)
        
        x=np.concatenate([[det_left_width, det_left_width, 0, 0, 0.5, 0.5, 0], -holeX-np.cos(phi)*holeR, [0, 0.5, 0.5, 0,0, det_left_width]])+gap/2
        y=np.concatenate([[-det_height/2,det_height/2, det_height/2,28.14, 28.14,18.14, 18.14],holeR*np.sin(phi),[-18.14, -18.14,-28.14, -28.14,-det_height/2, -det_height/2]])
        return Polygon(zip(x,y))
    
#thickness of the cell walls in the 3d-printed frame
wall_thickness = 0.08
#gap between the cell cell of the 3d printed frame and the scintillator tile.  
wall_scint_gap=0.01

nW=0
absorber_material=['W']*nW+['Fe']*(n_layers-nW+1)