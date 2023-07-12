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
def get_xyzr_reco_reweighted_H3(arrays, event, w0=4, weight_by_granularity=True, b=3, prefix="ZDC"):
    x=arrays[f'{prefix}HitsReco.position.x'][event]
    y=arrays[f'{prefix}HitsReco.position.y'][event]
    z=arrays[f'{prefix}HitsReco.position.z'][event]
    E=arrays[f'{prefix}HitsReco.energy'][event]
    #lay=arrays['HcalEndcapPInsertHitsReco.layer'][event]
    sl=arrays[f'{prefix}HitsReco.dimension.x'][event]/2
    
    tmp=sorted(list(set(z)))
    dz=tmp[1]-tmp[0]
    
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
        Eneighbors=[0,0,0,0,0,0]
        for j in range(len(x)):
            if abs(z[i]-z[j])>dz*1.1 or  E[j]<thresh or  j == i :
                continue
            dx=(x[j]-x[i])/sl[i]
            dy=(y[j]-y[i])/sl[i]
            if abs(dx)>1.1 or abs(dy)>1.1:
                continue
            tol=0.01
            for k in range(6):
                if abs(dx-cph[k])<tol and abs(dy-sph[k])<tol:
                    #print("found neighbor")
                    Eneighbors[k]=E[j]
                    #print(Eneighbors, E[j])
                    break
        a=thresh*b
        Eneighbors=np.array(Eneighbors)
        #print(Eneighbors)
        reweight_energy=mx(Eneighbors,a)*mx(np.roll(Eneighbors,1),a)
        reweight_energy/=sum(reweight_energy)
        #print(reweight_energy)
        for k in range(6):
            xnew.append(x[i]+sl[i]*(cph[k]+cph[k-5])/3)
            ynew.append(y[i]+sl[i]*(sph[k]+sph[k-5])/3)
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

if __name__ == "__main__":
    import sys
    infile=sys.argv[1]
    outfile=sys.argv[2]
    useH3reweighting=False
    if len(sys.argv)>3:
        if sys.argv[3] == "H3":
            useH3reweighting=True
    arrays=ur.open(f'{infile}:events').arrays()
    w0=5; b=0.05
    w0_nrw=4
    drs=[]
    drs_rw=[]
    dxs=[]
    dys=[]
    dxs_rw=[]
    dys_rw=[]
    Es=[]
    nevents=300
    for event in range(nevents):
        #print(arrays['ZDCHitsReco.energy'][event])
        #print(len(arrays['ZDCHitsReco.position.x'][event]))
        try:

            x_reco, y_reco, _, r_reco=get_xyzr_reco_no_reweighting(arrays, event, w0=w0_nrw, weight_by_granularity=True, prefix='ZDC')
            if useH3reweighting:
                x_reco_rw, y_reco_rw, _, r_reco_rw=get_xyzr_reco_reweighted_H3(arrays, event, w0=w0, b=b, weight_by_granularity=True, prefix='ZDC')
            x_truth, y_truth, _, r_truth=get_xyzr_truth(arrays, event, w0=w0, weight_by_granularity=True, prefix='ZDC')
            #print(r_truth, r_reco)
            drs.append(r_reco-r_truth)
            dxs.append(x_reco-x_truth)
            dys.append(y_reco-y_truth)
            if useH3reweighting:
                drs_rw.append(r_reco_rw-r_truth)
                dxs_rw.append(x_reco_rw-x_truth)
                dys_rw.append(y_reco_rw-y_truth)
            Es.append(sum(arrays[f'ZDCHitsReco.energy'][event]))
        except:
            pass

        if event%10==0:
            print(f"{infile}: done with event {event}/{nevents}")
    d=dict(E=Es, dr=drs, dy=dys, dx=dxs)
    if useH3reweighting:
        d["dr_rw"]=drs_rw
        d["dx_rw"]=dxs_rw
        d["dy_rw"]=dys_rw
    print({a:len(d[a]) for a in d})
    pd.DataFrame(d).to_csv(outfile)
