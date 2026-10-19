import scipy.io as sio, numpy as np
D='/tmp/claude-0/-home-user-Icode/1f385af9-cc34-5eb3-921b-4648ebd049a8/scratchpad/data/'
SETS={'CS':(D+'851347a7-CS_P_mat_Em_mat_Obs_mat.mat','CS_RealFinaldatacell','CSseq'),
      'fumin':(D+'e6c6b609-shortF__extradays_P_mat_Em_mat_Obs_data.mat','RealFinaldatacell','fuminobservable'),
      'sss':(D+'a7318b9f-shortS__P_mat_Em_mat_Obs_data.mat','sss_RealFinaldatacell','Sssseq'),
      'per0':(D+'360d9e0a-per0__P_mat_Em_mat_Obs_data.mat','per0_RealFinaldatacell','Per0seq')}
def windows(g):
    f,ck,rk=SETS[g]; m=sio.loadmat(f); C=m[ck]; R=m[rk]; out=[]
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            cm=C[3*i,j]; y=R[i,j]
            if np.size(cm)!=1 or np.size(y)==0 or int(cm.ravel()[0])!=3: continue
            P=np.asarray(C[3*i+1,j],float); B=np.asarray(C[3*i+2,j],float)
            if P.shape!=(3,3) or B.shape!=(3,2): continue
            out.append(dict(fly=i+1,win=j+1,P=P,B=B,y=np.asarray(y).ravel().astype(int)-1))
    return out
def xn(P):
    kX=P[...,0,1]/(P[...,0,1]+P[...,0,2]); a=P[...,2,0]*kX; return (a-P[...,2,1])/(a+P[...,2,1])
def yn(P):
    kY=P[...,0,2]/(P[...,0,1]+P[...,0,2]); a=P[...,1,0]*kY; return (a-P[...,1,2])/(a+P[...,1,2])
