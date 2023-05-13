from shapely import *

import beampipe_parameters
beampipe0=beampipe_parameters.Beampipe()
import numpy as np

def layer_boundaries(layer=0, side="L", beampipe=beampipe0, height=60.96, width=60.96,gap_between_sides=0.4,
                    gap_for_backplanes=0.4, backplanes_in_holes=True):
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