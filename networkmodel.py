'''
Created on 2019/11/17

@author: Takagi
'''
import numpy as np



def setLineProb(pi3, inp0):
    b=2
    for i in range(3):
        n=i+b-1
        inp=np.copy(inp0)
        inp+=b
        inp=np.where(inp!=n,0,inp)
        inp=np.where(inp==n,1,inp)            
        pi3[:,i]=np.sum(inp,axis=0)


def getMatrixAtval(inp0,i):
    b=2
    n=i+b-1
    inp=np.copy(inp0)
    inp+=b
    inp=np.where(inp!=n,0,inp)
    inp=np.where(inp==n,1,inp)            
    return inp

def prodProduct(inp,oup):
    mat=np.dot(inp.T,oup)
    return mat

def setMutProb(mu3,inp0,oup0):
    sinps=[]
    soups=[]
    for i in range(3):
        sinp=getMatrixAtval(inp0,i)
        soup=getMatrixAtval(oup0,i)
        sinps.append(sinp)
        soups.append(soup)
    for i in range(3):
        for j in range(3):
            k=i*3+j
            mu3[:,:,k]=prodProduct(sinps[i],soups[j])

def prep_state_probability(inp, oup):
    numcomp=0
    if inp.shape[0]>0:
        numcomp=inp.shape[1]
    pi3=np.zeros((numcomp,3))
    p3=np.zeros((numcomp,3))
    setLineProb(pi3, inp)
    setLineProb(p3, oup)       
    mu3=np.zeros((numcomp,numcomp,9))
    
    inp.astype(np.int)
    oup.astype(np.int)
    setMutProb(mu3,inp,oup)        
    return mu3, pi3, p3


def Normalize3StLine( pi3, p3):
    t1=np.zeros(pi3.shape[0])
    t2=np.zeros(p3.shape[0])
    for i in range(pi3.shape[1]):
        t1+=pi3[:,i]
        t2+=p3[:,i]
    t1=np.where(t1==0,1,t1)
    t2=np.where(t2==0,1,t2)
    for i in range(pi3.shape[1]):
        pi3[:,i]=pi3[:,i]/t1
        p3[:,i]=p3[:,i]/t2


def Normalize3StMat(mu3):
    t=np.zeros(mu3.shape[0]*mu3.shape[1])
    t=np.reshape(t,(mu3.shape[0], mu3.shape[1]))
    for z in range(mu3.shape[2]):
        t+=mu3[:,:,z]
    t=np.where(t==0,1,t)
    for z in range(mu3.shape[2]):
        mu3[:,:,z]=mu3[:,:,z]/t

def normaliZe3State( mu3, pi3, p3):
    Normalize3StLine(pi3, p3)
    Normalize3StMat(mu3)

def LineEntrop(pi3):
#        pi3=np.copy(pi30)
    pi3=np.where(pi3==0,1,pi3)
    s10=pi3*(-np.log(pi3))
    s1=np.sum(s10)
    return s1
    

def getAllMutJoiInfo3StatDef(mu30, pi30, p30):
    mu3=np.copy(mu30)
    pi3=np.copy(pi30)
    p3=np.copy(p30)
    x=mu3.shape[0]
    y=mu3.shape[1]
    s1=x*LineEntrop(pi3)
    s2=y*LineEntrop(p3)
    s3=LineEntrop(mu3)
    ret=s1+s2-s3
    retn=x
    if retn==0:
        retn=1
    ret= ret/retn
    l2= np.log(2.0)
    ret=ret/l2
    return ret
    
def estimate_InOut_signal_mutal_entrop(inp0, oup0):
    inp=np.copy(inp0)
    oup=np.copy(oup0)
    mu3, pi3, p3=prep_state_probability(inp, oup)
    normaliZe3State(mu3, pi3, p3)        
    ret=getAllMutJoiInfo3StatDef(mu3, pi3, p3)
    
    return ret


def get_node_activity_cost(inp, oup, mat):
    celln=oup.shape[1]
    cout=np.dot(inp, mat)
    cce=cout*oup
    ret=np.zeros(celln)
    for j in range(celln):
        ret[j]=np.sum(cce[:,j])
    
    tot=oup.shape[0]
    if tot==0:
        tot=1
    ret=ret/tot
    return ret 

def set_ouput_signals(pdat,sdv):
    th1=sdv
    th2=-sdv                
    pdat=np.where((pdat<th1) &(pdat>th2),0,pdat)
    pdat=np.where((pdat>0),1,pdat)
    pdat=np.where((pdat<0),-1,pdat)
    return pdat
 

def get_output_signal(mat, th, signal_input):
    prod=np.dot(signal_input,mat)
    psig= set_ouput_signals(prod, th)
    return psig
    


def get_InOut_signal(mat, th, signal_input):
    siz=mat.shape[0]
    output = get_output_signal(mat, th, signal_input)
    oup=np.array(output)
    return signal_input, oup


def set_matrix_threshold_pn(dat,fth):
    mat_th=np.where((dat>=-fth) & (dat<=fth),0,dat)
    return mat_th


def get_avr_sd(pd):
    av= np.average(pd, None, None, False)
    sd=np.std(pd)
    ret= np.zeros(2)
    ret[0]=av
    ret[1]=sd
    return ret

def set_matrix_thershold(dat0, th):
    dat=np.copy(dat0)
    datab= np.abs(dat)
    avsd= get_avr_sd(datab)
    thresh_val= avsd[0]+avsd[1]*th
    mat_th=set_matrix_threshold_pn(dat, thresh_val)
    return mat_th, thresh_val


def getOutputOverlap(pat, matth):
    tot=0
    res=0
    for i in range(len(pat)):
        for j in range(len(pat)):
            if i!=j:
                df=pat[i]-pat[j]
                df=np.where(df!=0,1,df)
                df=np.sum(df)
                res+=df
                tot+=1
    if tot==0:
        tot=1
    res=res/tot
    return res    

def get_network_energy(dt,th,signal_input):
        matth, thresh_val=set_matrix_thershold(dt,th)
        mat1=np.copy(matth)
        mat1=np.where(mat1!=0,1,mat1)        
        wiring=np.sum(matth!=0)
        #wiring=np.sum(mat1!=0)
        size=dt.shape[0]
        if size==0:
            size=1
        activity_list=[]
        overlap_list=[]
        mut_entropy_list=[]
        for sig in signal_input:
            inp0, oup0=get_InOut_signal(matth,thresh_val,sig)
            node_ene=get_node_activity_cost(np.abs(oup0), np.abs(oup0),np.abs(matth))
            activity_cost=np.sum(node_ene)
            activity_list.append(activity_cost)

            single_mut_ent=estimate_InOut_signal_mutal_entrop(inp0, oup0)
            mut_entropy_list.append(single_mut_ent)
            single_overlap=getOutputOverlap(oup0, matth)
            overlap_list.append(single_overlap)
            
        activity=np.average(activity_list)
        overlap=np.average(overlap_list)
        overlap=(size-overlap)/size
        mut_entropy=np.average(mut_entropy_list)
        
        if activity!=0:
            wiring_par_activity=wiring/activity
        else:
            wiring_par_activity=0
        ret=[mut_entropy,overlap,wiring_par_activity]
        return ret    

def get_simulation_result(thres, dt, signal_input):
    ret=[]
    for th in thres:
        pret=get_network_energy(dt,th,signal_input)
        ret.append(pret)
        
    return ret    


def read_connectivity_data(file):
    d= np.loadtxt(file)
    return d

def createRandomInitial(matrix_size, repeat_num, activation_probability):
    pdat=np.random.rand(repeat_num*matrix_size)
    pdat=np.array(pdat)
    pdat=np.reshape(pdat,(repeat_num,matrix_size))
    th1=1-activation_probability
    th2=activation_probability
    pdat=np.where((pdat<th1) &(pdat>th2),0,pdat)
    pdat=np.where((pdat>th1),1,pdat)
    pdat=np.where((pdat<th2) & (pdat>0),-1,pdat)
    return pdat

def create_random_input(batch_size,matrix_size,repeat_num, activation_probability):
    ret=[]
    for i in range(batch_size):
        sig=createRandomInitial(matrix_size, repeat_num, activation_probability)
        ret.append(sig)
    return ret


def checkPath(path):
    import os
    if not os.path.exists(path):
        os.mkdir(path)


def checkinputfolder(input_folder):
    import glob    
    files=glob.glob(input_folder+"/*.txt")
    #input connectivity file list (txt)
    return files
def getThresholds(sti, edi,threshold_step):
    #threshold_step=0.1
    n=0
    pp=sti
    ret=[]
    while pp<=edi:
        pp=n*threshold_step+sti
        ret.append(pp)
        n+=1
    return ret


def print_out_result_test(ret):
    print("mut_entropy \t overlap \t wiring_par_activity")
    for r in ret:
        s=""
        for v in r:
            s+="{}\t".format(v)
        print(s)

def networksimulation(set):
    #execute simulation
    #set=setting.py
    files=checkinputfolder(set.input_folder)
    checkPath(set.output_folder)
    if len(files)>0:
        thrs=getThresholds(set.start_threshold,set.end_threshold,set.threshold_step )
        signal_input=create_random_input(set.batch_size,set.matrix_size,set.repeat_num, set.activation_probability)
        for file in files:
            if set.deb==1:
                print("Input file:"+file)
            dt=read_connectivity_data(file)
            ret=get_simulation_result(thrs, dt, signal_input)
            if set.deb==1:
                print_out_result_test(ret)






