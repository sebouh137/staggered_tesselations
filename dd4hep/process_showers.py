import uproot as ur, numpy as np, pandas as pd

c_in_mm_per_ns = 299.792458

#default values.  Can be configured
Emin=0.000472*0.1 # 1/10 a MIP
tmax=200 # 200 ns

def get_xyzr_reco_no_reweighting(arrays, event, w0=4, weight_by_granularity=True, prefix="ZDC", MIP=0.00047):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    t=arrays[f'{prefix}HitsReco.time'][event] - z/c_in_mm_per_ns #correct for time of flight
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    #make selection of hits to be used
    slc=(E>Emin) & (t<tmax)
    x=x[slc]
    y=y[slc]
    z=z[slc]
    E=E[slc]
    sl=sl[slc]

    if type(w0)!=float or w0 !=0 :
        w=w0+np.log((E+.0000001)/sum(E))
        w=w*(w>0)
    else :
        w=E
    if weight_by_granularity:
        w=w/sl**2
    
    sumw=np.sum(w+.0000001, axis=-1)
    x_reco=np.sum(x*w, axis=-1)/sumw
    y_reco=np.sum(y*w, axis=-1)/sumw
    z_reco=np.sum(z*w, axis=-1)/sumw
    r_reco=np.hypot(x_reco,y_reco)
    return [x_reco,y_reco,z_reco,r_reco]

def get_xyzr_truth(arrays, event, w0=4, weight_by_granularity=True, prefix="ZDC"):
    #first determine the z value used in the recon
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    t=arrays[f'{prefix}HitsReco.time'][event] - z/c_in_mm_per_ns #correct for time of flight
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]
    
    #make selection
    slc=(E>Emin) & (t<tmax)
    z=z[slc]
    E=E[slc]
    t=t[slc]
    sl=sl[slc]
    
    if type(w0)!= float or w0 !=0:
        w=w0+np.log((E+.0000001)/sum(E))
        w=w*(w>0)
    else :
        w=E
    if weight_by_granularity:
        w=w/sl**2
    
    z_reco=np.sum(z*w, axis=-1)/np.sum(w+.000000001, axis=-1)
    
    px=arrays["MCParticles.momentum.x"][event,2]
    py=arrays["MCParticles.momentum.y"][event,2]
    pz=arrays["MCParticles.momentum.z"][event,2]

    x_truth=px/pz*z_reco
    y_truth=py/pz*z_reco
    z_truth=z_reco
    r_truth=np.hypot(x_truth,y_truth)
    
    return [x_truth,y_truth,z_truth,r_truth]

sqrt3=np.sqrt(3)
sqrt5=np.sqrt(5)

def mx(a,b):
    return a*(a>b)+b*(b>=a)
#with H3
def get_xyzr_reco_reweighted_H3(arrays, event, w0=5, weight_by_granularity=True, prefix="ZDC", MIP=0.000472):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event] 
    E=arrays[f'{prefix}HitsReco.energy'][event]
    #lay=arrays['HcalEndcapPInsertHitsReco.layer'][event]
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    t=arrays[f'{prefix}HitsReco.time'][event]  - z/c_in_mm_per_ns #correct for time of flight
    slc=(E>Emin) & (t<tmax)
    x=x[slc]
    y=y[slc]
    z=z[slc]
    E=E[slc]
    sl=sl[slc]

    
    minz=min(z)
    dz=min(z[z!=minz])-minz
    
    Etot=sum(E)
    if type(w0) !=float or  w0 != 0.0:
        thresh=Etot*np.exp(-np.max(w0))
    else :
        thresh = 0
    #print("thresh", thresh)
    
    xnew=[]
    ynew=[]
    znew=[]
    Enew=[]
    slnew=[]
    phi=np.linspace(0, np.pi*5/3, 6)
    cph=np.cos(phi)
    sph=np.sin(phi)
    roll_cph=np.roll(cph,1)
    roll_sph=np.roll(sph,1)
    
    #centers of the triangles minus that of the hexagon
    dx_tri=(cph+roll_cph)/3
    dy_tri=(sph+roll_sph)/3
    
    for i in range(len(x)):
        if E[i]<thresh:
            continue
        neighbors_found=0
        Eneighbors=[0,0,0,0,0,0]
        for j in range(len(x)):

            if abs(z[i]-z[j])>dz*1.1 or E[j]<thresh or z[i]==z[j]:
                continue
            dx=(x[j]-x[i])/sl[i]
            dy=(y[j]-y[i])/sl[i]
            if abs(dx)>1.1 or abs(dy)>1.1:
                continue
            tol=0.01

            
            if neighbors_found==6:
                break
            
            for k in range(6):
                if Eneighbors[k]:
                    continue
                if abs(dx-cph[k])<tol and abs(dy-sph[k])<tol:
                    #print("found neighbor")
                    Eneighbors[k]=E[j]
                    #print(Eneighbors, E[j])
                    neighbors_found+=1
                    break
        #a=thresh*b
        
        Eneighbors=np.array(Eneighbors)

        reweight_energy=mx(Eneighbors,MIP)*mx(np.roll(Eneighbors,1),MIP) 
        reweight_energy/=sum(reweight_energy)

        for k in range(6):
            if E[i]*reweight_energy[k]<thresh:
                continue
            xnew.append(x[i]+sl[i]*dx_tri[k])
            ynew.append(y[i]+sl[i]*dy_tri[k])
            znew.append(z[i])
            Enew.append(E[i]*reweight_energy[k])
            slnew.append(sl[i]/sqrt5)
            
    xnew=np.array(xnew)
    ynew=np.array(ynew)
    znew=np.array(znew)
    Enew=np.array(Enew)
    slnew=np.array(slnew)
    
    
    if type(w0)!=float or  w0 !=0:
        w=w0+np.log((Enew+.0000001)/Etot)
        w=w*(w>0)
    else :
        w=Enew
    if weight_by_granularity:
        w=w/slnew**2

    sumw=np.sum(w+.0000001, axis=-1)
    x_reco=np.sum(xnew*w, axis=-1)/sumw
    y_reco=np.sum(ynew*w, axis=-1)/sumw
    z_reco=np.sum(znew*w, axis=-1)/sumw
    r_reco=np.hypot(x_reco,y_reco)
    return [x_reco,y_reco,z_reco,r_reco]

#with H4                                                                                                     
def get_xyzr_reco_reweighted_H4(arrays, event, w0=6, weight_by_granularity=True, prefix="ZDC", MIP=0.000472):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    t=arrays[f'{prefix}HitsReco.time'][event] - z/c_in_mm_per_ns #correct for time of flight
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    # time and energy cuts
    slc=(E>Emin) & (t<tmax)
    x=x[slc]
    y=y[slc]
    z=z[slc]
    E=E[slc]
    sl=sl[slc]
    

    minz=min(z)
    dz=min(z[z!=minz])-minz

    Etot=sum(E)
    if type(w0) !=float or  w0 != 0.0:
        thresh=Etot*np.exp(-np.max(w0))
    else :
        thresh = 0


    xnew=[]
    ynew=[]
    znew=[]
    Enew=[]
    slnew=[]
    phi=np.linspace(0, np.pi*5/3, 6)
    cph=np.cos(phi)
    sph=np.sin(phi)

    for i in range(len(x)):
        if E[i]<thresh:
            continue
        # there are twelve neighboring cell positions where we need to determine the energy
        Eneighbors=[0,0,0,0,0,0,0,0,0,0,0,0]
        for j in range(len(x)):
            if abs(z[i]-z[j])>dz*2.1 or  E[j]<thresh or  j == i :
                continue
            dx=(x[j]-x[i])/sl[i]
            dy=(y[j]-y[i])/sl[i]
            if abs(dx)>1.6 or abs(dy)>1.6:
                continue
            tol=0.01
            for k in range(6):
                if abs(dx-1.5*cph[k])<tol and abs(dy-1.5*sph[k])<tol:
                    Eneighbors[k]+=E[j]
                    break
            for k in range(6):
                if abs(dx+sqrt3/2*sph[k])<tol and abs(dy-sqrt3/2*cph[k])<tol:
                    Eneighbors[k+6]+=E[j]
                    break
        
        Eneighbors=mx(np.array(Eneighbors),MIP)
        
        
        reweight_energy_1=Eneighbors[:6]*np.roll(Eneighbors[6:],4)*np.roll(Eneighbors[6:],5)
        reweight_energy_2=Eneighbors[6:]*np.roll(Eneighbors[6:],-1)*np.roll(Eneighbors[6:],1)
        reweight_energy=np.concatenate([reweight_energy_1, reweight_energy_2])
        reweight_energy/=sum(reweight_energy)

        for k in range(6):
            if E[i]*reweight_energy[k]< thresh:
                continue
            xnew.append(x[i]+sl[i]*0.75*cph[k])
            ynew.append(y[i]+sl[i]*0.75*sph[k])
            znew.append(z[i])
            Enew.append(E[i]*reweight_energy[k])
            slnew.append(sl[i]/sqrt5)
        
        for k in range(6):
            if E[i]*reweight_energy[k+6]<thresh:
                continue
            xnew.append(x[i]+sl[i]*sqrt3/4*-sph[k])
            ynew.append(y[i]+sl[i]*sqrt3/4*cph[k])
            znew.append(z[i])
            Enew.append(E[i]*reweight_energy[k+6])
            slnew.append(sl[i]/sqrt5)
        
    
    xnew=np.array(xnew)
    ynew=np.array(ynew)
    znew=np.array(znew)
    Enew=np.array(Enew)
    
    slnew=np.array(slnew)
        
    if type(w0)!=float or  w0 !=0:
        w=w0+np.log((Enew+.0000001)/Etot)
        w=w*(w>0)
    else :
        w=Enew
    if weight_by_granularity:
        w=w/slnew**2

    sumw=np.sum(w+.0000001, axis=-1)
    x_reco=np.sum(xnew*w, axis=-1)/sumw
    y_reco=np.sum(ynew*w, axis=-1)/sumw
    z_reco=np.sum(znew*w, axis=-1)/sumw
    r_reco=np.hypot(x_reco,y_reco)

    return [x_reco,y_reco,z_reco,r_reco]

sqrt2=np.sqrt(2)
#with S2
def get_xyzr_reco_reweighted_S2(arrays, event, w0=6, weight_by_granularity=True, prefix="ZDC", MIP=0.000472):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event] - z/c_in_mm_per_ns #correct for time of flight
    E=arrays[f'{prefix}HitsReco.energy'][event]
    
    #side length.  
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    t=arrays[f'{prefix}HitsReco.time'][event]
    slc=(E>Emin) & (t<tmax)
    x=x[slc]
    y=y[slc]
    z=z[slc]
    E=E[slc]
    sl=sl[slc]

    
    minz=min(z)
    dz=min(z[z!=minz])-minz

    Etot=sum(E)
    if type(w0) !=float or  w0 != 0.0:
        thresh=Etot*np.exp(-np.max(w0))
    else :
        thresh = 0


    xnew=[]
    ynew=[]
    znew=[]
    Enew=[]
    slnew=[]
    phi=np.linspace(np.pi/4, np.pi*7/4, 4)
    cph=np.cos(phi)
    sph=np.sin(phi)
    #print("bbb")
    for i in range(len(x)):
        if E[i]<thresh:
            continue
        neighbors_found = 0
        Eneighbors=[0,0,0,0]
        for j in range(len(x)):
                                
            if abs(z[i]-z[j])>dz*1.1 or  E[j]<thresh or  j == i or z[i]==z[j] :
                continue
            dx=(x[j]-x[i])/sl[i]
            dy=(y[j]-y[i])/sl[i]
            if abs(dx)>0.6 or abs(dy)>0.6:
                continue
            tol=0.01
            for k in range(4):
                if abs(dx-sqrt2*cph[k]/2)<tol and abs(dy-sqrt2*sph[k]/2)<tol:
                    Eneighbors[k]+=E[j]
                    neighbors_found+=1
                    break
            if neighbors_found == 8:
                break
            
        #print("??")

        Eneighbors=mx(np.array(Eneighbors),MIP)
        #print("Eneighbors:", Eneighbors)
        
        reweight_energy=Eneighbors/sum(Eneighbors)
        #print("reweight_energy:", reweight_energy)
        #print(reweight_energy)
        #print(Eneighbors)
        #reweight_energy/=sum(reweight_energy)

        #print("aaa")
        for k in range(4):
            if E[i]*reweight_energy[k]< thresh:
                continue
            xnew.append(x[i]+sl[i]*sqrt2/4*cph[k])
            ynew.append(y[i]+sl[i]*sqrt2/4*sph[k])
            znew.append(z[i])
            Enew.append(E[i]*reweight_energy[k])
            slnew.append(sl[i]/2)
        
    
    xnew=np.array(xnew)
    ynew=np.array(ynew)
    znew=np.array(znew)
    Enew=np.array(Enew)
    
    slnew=np.array(slnew)
    #print("fff")
        
    if type(w0)!=float or  w0 !=0:
        w=w0+np.log((Enew+.0000001)/Etot)
        w=w*(w>0)
    else :
        w=Enew
    if weight_by_granularity:
        w=w/slnew**2
    #print("ggg")
    sumw=np.sum(w+.0000001, axis=-1)
    x_reco=np.sum(xnew*w, axis=-1)/sumw
    y_reco=np.sum(ynew*w, axis=-1)/sumw
    z_reco=np.sum(znew*w, axis=-1)/sumw
    r_reco=np.hypot(x_reco,y_reco)
    #print("hh")
    return [x_reco,y_reco,z_reco,r_reco]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(
                    prog='process_showers',
                    description='determines the position of the showers',
                    epilog='Text at the bottom of help')
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('--H4', help="flag for using H4 HEX-SPLIT",
                    action='store_true')  # on/off flag
    parser.add_argument('--H3',	help="flag for using H3 HEX-SPLIT",
                    action='store_true')  # on/off flag 
    parser.add_argument('--S2',  help="flag for using S2 HEX-SPLIT",
                    action='store_true')  # on/off flag
    parser.add_argument('-n', '--nevents', help="number of events to run", default=-1, type=int)
    parser.add_argument("-p", '--prefix', help="prefix for the detector type (hit type)", default="ZDC")
    parser.add_argument("-s", '--skip', help="number of events to skip", default=0, type=int)
    parser.add_argument("--w0_rw", help="w0 used in reweighted", default=5, type=float)
    parser.add_argument("--w0_nrw", help="w0 used in reweighted", default=4, type=float)
    parser.add_argument("--w0_use_range", help="use a range for w0", action="store_true")
    parser.add_argument("--MIP", help="MIP value in GeV", default=0.000472, type=float)
    args = parser.parse_args()
    
    import sys
    infile=args.infile
    outfile=args.outfile
    useH3reweighting=args.H3
    useH4reweighting=args.H4
    useS2reweighting=args.S2
    nevents=args.nevents
    prefix=args.prefix
    first_event=args.skip
    arrays=ur.open(f'{infile}:events').arrays()
    w0=args.w0_rw
    MIP=args.MIP
    w0_nrw=args.w0_nrw

    w0_use_range=args.w0_use_range
    if w0_use_range:
        w0_nrw=np.array([[a] for a in np.linspace(3.0, 8.0, 21)])
        w0=np.array([[a] for a in np.linspace(3.0, 8.0, 21)])
    x_truths=[]
    y_truths=[]
    drs=[]
    drs_rw=[]
    dxs=[]
    dys=[]
    dxs_rw=[]
    dys_rw=[]
    Es=[]
    mc_pzs=[]
    w0s=[]
    if nevents==-1:
        nevents=len(arrays)
    
    for event in range(first_event,first_event+nevents):
        #print(arrays['ZDCHitsReco.energy'][event])
        #print(len(arrays['ZDCHitsReco.position.x'][event]))
        try:

            x_reco, y_reco, _, r_reco=get_xyzr_reco_no_reweighting(arrays, event, w0=w0_nrw, weight_by_granularity=True, prefix=prefix)
            if useH3reweighting:
                x_reco_rw, y_reco_rw, _, r_reco_rw=get_xyzr_reco_reweighted_H3(arrays, event, w0=w0, MIP=MIP, weight_by_granularity=True, prefix=prefix)
            elif useH4reweighting:
                x_reco_rw, y_reco_rw, _, r_reco_rw=get_xyzr_reco_reweighted_H4(arrays, event, w0=w0, MIP=MIP, weight_by_granularity=True, prefix=prefix)
            elif useS2reweighting:
                x_reco_rw, y_reco_rw, _, r_reco_rw=get_xyzr_reco_reweighted_S2(arrays, event, w0=w0, MIP=MIP, weight_by_granularity=True, prefix=prefix)

            x_truth, y_truth, _, r_truth=get_xyzr_truth(arrays, event, w0=w0, weight_by_granularity=True, prefix=prefix)
            #print(r_truth, r_reco)
            drs.append(r_reco-r_truth)
            dxs.append(x_reco-x_truth)
            dys.append(y_reco-y_truth)
            x_truths.append(x_truth)
            y_truths.append(y_truth)
            if useH3reweighting or useH4reweighting or useS2reweighting:
                drs_rw.append(r_reco_rw-r_truth)
                dxs_rw.append(x_reco_rw-x_truth)
                dys_rw.append(y_reco_rw-y_truth)
            Es.append(sum(arrays[f'{prefix}HitsReco.energy'][event]))
            mc_pzs.append(arrays[f'MCParticles.momentum.z'][event,2])
            if w0_use_range:
                w0s.append(w0)
        except Exception as e:
            print(e)
            pass

        if event%10==0:
            print(f"{infile}: done with event {event}/{nevents}")
    if not w0_use_range:
        d=dict(E=Es, dr=drs, dy=dys, dx=dxs, mc_pz=mc_pzs, x_truth=x_truths, y_truth=y_truths)
        if useH3reweighting or useH4reweighting or useS2reweighting:
            d["dr_rw"]=drs_rw
            d["dx_rw"]=dxs_rw
            d["dy_rw"]=dys_rw
    else :
        w0s = [a[0] for a in w0] # flatten array
        d=dict(E=Es, dr=drs, dy=dys, dx=dxs, mc_pz=mc_pzs, x_truth=x_truths, y_truth=y_truths)
        if useH3reweighting or useH4reweighting or useS2reweighting:
            d["dr_rw"]=drs_rw
            d["dx_rw"]=dxs_rw
            d["dy_rw"]=dys_rw
        for key in list(d.keys()):
            #print(key)
            #print(type(d[key]))
            #print(d[key][0])
            if key not in "E mc_pz w0s".split():
                print("making new columns")
                for i in range(len(w0s)):
                    d[f"{key}_w0_{w0s[i]}".replace(".", "pt")]=[d[key][j][i] for j in range(len(d[key]))]
            
    
    print({a:len(d[a]) for a in d})
    if ".csv" in outfile:
        pd.DataFrame(d).to_csv(outfile)
    if ".pkl" in outfile:
        pd.DataFrame(d).to_pickle(outfile)
