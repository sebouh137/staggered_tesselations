import uproot as ur, numpy as np, pandas as pd

def get_xyzr_reco_no_reweighting(arrays, event, w0=4, weight_by_granularity=True, prefix="ZDC"):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    #lay=arrays['HcalEndcapPInsertHitsReco.layer'][event]
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    if w0 !=0:
        w=w0+np.log((E+.0000001)/sum(E))
        w=w*(w>0)
    else :
        w=E
    if weight_by_granularity:
        w=w/sl**2
    sumw=sum(w)
    x_reco=sum(x*w)/sumw
    y_reco=sum(y*w)/sumw
    z_reco=sum(z*w)/sumw
    r_reco=np.hypot(x_reco,y_reco)
    return [x_reco,y_reco,z_reco,r_reco]

def get_xyzr_truth(arrays, event, w0=4, weight_by_granularity=True, prefix="ZDC"):
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    #lay=arrays['HcalEndcapPInsertHitsReco.layer'][event]
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]
    if w0 !=0:
        w=w0+np.log((E+.0000001)/sum(E))
        w=w*(w>0)
    else :
        w=E
    if weight_by_granularity:
        w=w/sl**2
    z_reco=sum(z*w)/sum(w)
    
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
def get_xyzr_reco_reweighted_H3(arrays, event, w0=5, weight_by_granularity=True, b=0.00001, prefix="ZDC"):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    #lay=arrays['HcalEndcapPInsertHitsReco.layer'][event]
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    
    minz=min(z)
    dz=min(z[z!=minz])-minz
    
    Etot=sum(E)
    thresh=Etot*np.exp(-w0)
    
    
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
        a=thresh*b
        Eneighbors=np.array(Eneighbors)

        reweight_energy=mx(Eneighbors,a)*mx(np.roll(Eneighbors,1),a) 
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
    
    
    if w0 !=0:
        w=w0+np.log((Enew+.0000001)/Etot)
        w=w*(w>0)
    else :
        w=Enew
    if weight_by_granularity:
        w=w/slnew**2
    #print(w, Enew/Etot/thresh)
    sumw=np.sum(w)
    x_reco=np.sum(xnew*w)/sumw
    y_reco=np.sum(ynew*w)/sumw
    z_reco=np.sum(znew*w)/sumw
    r_reco=np.hypot(x_reco,y_reco)
    return [x_reco,y_reco,z_reco,r_reco]

#with H4                                                                                                     
def get_xyzr_reco_reweighted_H4(arrays, event, w0=6, weight_by_granularity=True, b=0.0001, prefix="ZDC"):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]

    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2

    minz=min(z)
    dz=min(z[z!=minz])-minz

    Etot=sum(E)
    thresh=Etot*np.exp(-w0)


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
        Eneighbors=[0,0,0,0,0,0,0,0,0,0,0,0]
        for j in range(len(x)):
            if abs(z[i]-z[j])>dz*1.1 or  E[j]<thresh or  j == i :
                continue
            dx=(x[j]-x[i])/sl[i]
            dy=(y[j]-y[i])/sl[i]
            if abs(dx)>1.6 or abs(dy)>1.6:
                continue
            tol=0.01
            for k in range(6):
                if abs(dx-1.5*cph[k])<tol and abs(dy-1.5*sph[k])<tol:
                    Eneighbors[k]=E[j]
                    break
            for k in range(6):
                if abs(dx+sqrt3/2*sph[k])<tol and abs(dy-sqrt3/2*cph[k])<tol:
                    Eneighbors[k+6]=E[j]
                    break
        #print("??")
        a=thresh*b
        Eneighbors=mx(np.array(Eneighbors),a)
        #print("????")
        
        reweight_energy_1=Eneighbors[:6]*np.roll(Eneighbors[6:],4)*np.roll(Eneighbors[6:],5)
        reweight_energy_2=Eneighbors[6:]*np.roll(Eneighbors[6:],-1)*np.roll(Eneighbors[6:],1)
        reweight_energy=np.concatenate([reweight_energy_1, reweight_energy_2])
        #print(reweight_energy)
        #print(Eneighbors)
        reweight_energy/=sum(reweight_energy)


        for k in range(6):
            if E[i]*reweight_energy[k]< thresh:
                continue
            xnew.append(x[i]+sl[i]*0.75*cph[k])
            ynew.append(y[i]+sl[i]*0.75*sph[k])
            znew.append(z[i])
            Enew.append(E[i]*reweight_energy[k])
            slnew.append(sl[i]/sqrt5)
        #print("d")
        for k in range(6):
            if E[i]*reweight_energy[k+6]<thresh:
                continue
            xnew.append(x[i]+sl[i]*sqrt3/4*-sph[k])
            ynew.append(y[i]+sl[i]*sqrt3/4*cph[k])
            znew.append(z[i])
            Enew.append(E[i]*reweight_energy[k+6])
            slnew.append(sl[i]/sqrt5)
        #print("e")
    
    xnew=np.array(xnew)
    ynew=np.array(ynew)
    znew=np.array(znew)
    Enew=np.array(Enew)
    
    slnew=np.array(slnew)

    if w0 !=0:
        w=w0+np.log((Enew+.0000001)/Etot)
        w=w*(w>0)
    else :
        w=Enew
    if weight_by_granularity:
        w=w/slnew**2
    
    sumw=np.sum(w)
    x_reco=np.sum(xnew*w)/sumw
    y_reco=np.sum(ynew*w)/sumw
    z_reco=np.sum(znew*w)/sumw
    r_reco=np.hypot(x_reco,y_reco)
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

    parser.add_argument('-n', '--nevents', help="number of events to run", default=-1, type=int)
    parser.add_argument("-p", '--prefix', help="prefix for the detector type (hit type)", default="ZDC")
    parser.add_argument("-s", '--skip', help="number of events to skip", default=0, type=int)
    parser.add_argument("--w0_rw", help="w0 used in reweighted", default=5, type=float)
    args = parser.parse_args()
    
    import sys
    infile=args.infile
    outfile=args.outfile
    useH3reweighting=args.H3
    useH4reweighting=args.H4
    nevents=args.nevents
    prefix=args.prefix
    first_event=args.skip
    arrays=ur.open(f'{infile}:events').arrays()
    w0=args.w0_rw
    b=0.05
    w0_nrw=4

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

    if nevents==-1:
        nevents=len(arrays)
    
    for event in range(first_event,first_event+nevents):
        #print(arrays['ZDCHitsReco.energy'][event])
        #print(len(arrays['ZDCHitsReco.position.x'][event]))
        try:

            x_reco, y_reco, _, r_reco=get_xyzr_reco_no_reweighting(arrays, event, w0=w0_nrw, weight_by_granularity=True, prefix=prefix)
            if useH3reweighting:
                x_reco_rw, y_reco_rw, _, r_reco_rw=get_xyzr_reco_reweighted_H3(arrays, event, w0=w0, b=b, weight_by_granularity=True, prefix=prefix)
            elif useH4reweighting:
                x_reco_rw, y_reco_rw, _, r_reco_rw=get_xyzr_reco_reweighted_H4(arrays, event, w0=w0, b=b, weight_by_granularity=True, prefix=prefix)
            x_truth, y_truth, _, r_truth=get_xyzr_truth(arrays, event, w0=w0, weight_by_granularity=True, prefix=prefix)
            #print(r_truth, r_reco)
            drs.append(r_reco-r_truth)
            dxs.append(x_reco-x_truth)
            dys.append(y_reco-y_truth)
            x_truths.append(x_truth)
            y_truths.append(y_truth)
            if useH3reweighting or useH4reweighting:
                drs_rw.append(r_reco_rw-r_truth)
                dxs_rw.append(x_reco_rw-x_truth)
                dys_rw.append(y_reco_rw-y_truth)
            Es.append(sum(arrays[f'{prefix}HitsReco.energy'][event]))
            mc_pzs.append(arrays[f'MCParticles.momentum.z'][event,2])
        except:
            pass

        if event%10==0:
            print(f"{infile}: done with event {event}/{nevents}")
    d=dict(E=Es, dr=drs, dy=dys, dx=dxs, mc_pz=mc_pzs, x_truth=x_truths, y_truth=y_truths)
    if useH3reweighting or useH4reweighting:
        d["dr_rw"]=drs_rw
        d["dx_rw"]=dxs_rw
        d["dy_rw"]=dys_rw
    print({a:len(d[a]) for a in d})
    pd.DataFrame(d).to_csv(outfile)
