#!/usr/bin/python
# derive the minimum perturbation of CMAP left helix region to lower the population
# based on my fitcmap.py
# a large amount of rewrite to boost the speed
# the major principle is to minimize the calculation efforts in ecmap
# to seperate coor into one part which is related to observables (in this case simple)
# and another part which is related to ecmap (precompute as much as possible)
# March 2015, Jing Huang

import sys, string, copy
# import numpy as np # remove np so we can use pypy to speed up
import math
from math import sqrt, pi, cos
import random

def ecmapfast(ggrd, IP1, IP2, IP1P1, IP2P1, TT, TU):
# fast calculation of cmap energy 
    TX=[ggrd[0][IP2][IP1],ggrd[0][IP2][IP1P1],ggrd[0][IP2P1][IP1P1],ggrd[0][IP2P1][IP1], \
        ggrd[1][IP2][IP1],ggrd[1][IP2][IP1P1],ggrd[1][IP2P1][IP1P1],ggrd[1][IP2P1][IP1], \
        ggrd[2][IP2][IP1],ggrd[2][IP2][IP1P1],ggrd[2][IP2P1][IP1P1],ggrd[2][IP2P1][IP1], \
        ggrd[3][IP2][IP1],ggrd[3][IP2][IP1P1],ggrd[3][IP2P1][IP1P1],ggrd[3][IP2P1][IP1] ]

    inn=0
    TC=[[],[],[],[]]
    for i in range(4):
	for j in range(4):
            xx=0.0
            a=wt[inn]
	    for k in range(16):
		xx = xx + a[k]*TX[k]
	    inn+=1
	    TC[i].append(xx)

    E=0.0
    for II in [3,2,1,0]:
        E=TT*E+((TC[II][3]*TU+TC[II][2])*TU+TC[II][1])*TU+TC[II][0]
    return E

def convertvar(phi1,phi2):
    MRES=360.0/numspace
    iphi1=int(math.floor((phi1+180.0)/MRES))
    iphi2=int(math.floor((phi2+180.0)/MRES))

    TT=(phi1-iphi1*MRES+180.0)/MRES
    TU=(phi2-iphi2*MRES+180.0)/MRES

    IP1 = iphi1 % numspace  # charmm sometimes gives 180.0 dihedral angle
    IP2 = iphi2 % numspace
    IP1P1 = (IP1 + 1) % numspace
    IP2P1 = (IP2 + 1) % numspace
    return [IP1,IP2,IP1P1,IP2P1,TT,TU]

def ecmap(ggrd,phi1,phi2):
    MRES=360.0/numspace
    iphi1=int(math.floor((phi1+180.0)/MRES))
    iphi2=int(math.floor((phi2+180.0)/MRES))

    TY,TY1,TY2,TY12=gcstup2(ggrd,iphi1,iphi2)
    TC=gcscf(TY,TY1,TY2,TY12,MRES,MRES)
    TT=(phi1-iphi1*MRES+180.0)/MRES
    TU=(phi2-iphi2*MRES+180.0)/MRES

    E=0.0
    for II in [3,2,1,0]:
        E=TT*E+((TC[II][3]*TU+TC[II][2])*TU+TC[II][1])*TU+TC[II][0]
    return E

def readcmap(cmapfile):
    gmap=[]
    for i in range(numspace):
        gmap.append([])
    for i in range(numspace):
        cmapfile.readline()
        jcount=0
        for j in range(5):
            for str2 in string.split(cmapfile.readline()):        
                gmap[jcount].append(float(str2))
                jcount+=1
        cmapfile.readline()
    return gmap

def setcmap(cmap):
    # global ggrd
    ggrd=[]
    for i in range(4):
        ggrd.append(copy.deepcopy(cmap))  #ggrd[0], and create ggrd[1],ggrd[2],ggrd[3]
    num=numspace
    xm=num/2
    dnum=2*num # double of num, equal num+xm+xm
    dx=360.0/num
    xmin=-360.0

    tgmap=[]
    #expand i and j 
    for i in range(dnum):
        ii=(i+xm) % (num)
        tgmap.append([])
        for j in range(dnum):
            jj=(j+xm) % (num)
            tgmap[i].append(ggrd[0][jj][ii])
    #print tgmap

    y2a=[]
    for j in range(dnum):
        y2tmp=CMAPSPL(dx,tgmap[j],dnum)
        y2a.append(y2tmp)

    for ii in range(xm,num+xm):
        phi=(ii-xm)*dx-180.0
        for jj in range(xm,num+xm):
            psi=(jj-xm)*dx-180

            yytmp=[]
            y1tmp=[]
            for j in range(dnum):
                a,b=CMAPSPI(-360.0,dx,tgmap[j],y2a[j],psi)
                yytmp.append(a)
                y1tmp.append(b)
            # compute v1,v2,v12 for grid point phi,psi

            u2=CMAPSPL(dx,yytmp,dnum)
            v,v1=CMAPSPI(-360.0,dx,yytmp,u2,phi)
            u2=CMAPSPL(dx,y1tmp,dnum)
            v2,v12=CMAPSPI(-360.0,dx,y1tmp,u2,phi)
            ggrd[1][jj-xm][ii-xm]=v1
            ggrd[2][jj-xm][ii-xm]=v2
            ggrd[3][jj-xm][ii-xm]=v12
    return ggrd


def setcmap2(cmap):
    # notice that ggrd[1],ggrd[2],ggrd[3] is already multiplied by 15, 15, and 15*15
    # to be used with ecmapfast 
    # global ggrd
    ggrd=[]
    for i in range(4):
        ggrd.append(copy.deepcopy(cmap))  #ggrd[0], and create ggrd[1],ggrd[2],ggrd[3]
    num=numspace
    xm=num/2
    dnum=2*num # double of num, equal num+xm+xm
    dx=360.0/num
    xmin=-360.0

    tgmap=[]
    #expand i and j 
    for i in range(dnum):
        ii=(i+xm) % (num)
        tgmap.append([])
        for j in range(dnum):
            jj=(j+xm) % (num)
            tgmap[i].append(ggrd[0][jj][ii])
    #print tgmap

    y2a=[]
    for j in range(dnum):
        y2tmp=CMAPSPL(dx,tgmap[j],dnum)
        y2a.append(y2tmp)

    for ii in range(xm,num+xm):
        phi=(ii-xm)*dx-180.0
        for jj in range(xm,num+xm):
            psi=(jj-xm)*dx-180

            yytmp=[]
            y1tmp=[]
            for j in range(dnum):
                a,b=CMAPSPI(-360.0,dx,tgmap[j],y2a[j],psi)
                yytmp.append(a)
                y1tmp.append(b)
            # compute v1,v2,v12 for grid point phi,psi

            u2=CMAPSPL(dx,yytmp,dnum)
            v,v1=CMAPSPI(-360.0,dx,yytmp,u2,phi)
            u2=CMAPSPL(dx,y1tmp,dnum)
            v2,v12=CMAPSPI(-360.0,dx,y1tmp,u2,phi)
            ggrd[1][jj-xm][ii-xm]=v1*dx
            ggrd[2][jj-xm][ii-xm]=v2*dx
            ggrd[3][jj-xm][ii-xm]=v12*dx*dx
    return ggrd



def gcstup2(ggrd,P1,P2):
    # here setup evaluate the energy
    IP1 = P1 % numspace  # charmm sometimes gives 180.0 dihedral angle
    IP2 = P2 % numspace
    IP1P1 = (IP1 + 1) % numspace
    IP2P1 = (IP2 + 1) % numspace

    y=[ggrd[0][IP2][IP1],ggrd[0][IP2][IP1P1],ggrd[0][IP2P1][IP1P1],ggrd[0][IP2P1][IP1]]
    y1=[ggrd[1][IP2][IP1],ggrd[1][IP2][IP1P1],ggrd[1][IP2P1][IP1P1],ggrd[1][IP2P1][IP1]]
    y2=[ggrd[2][IP2][IP1],ggrd[2][IP2][IP1P1],ggrd[2][IP2P1][IP1P1],ggrd[2][IP2P1][IP1]]
    y12=[ggrd[3][IP2][IP1],ggrd[3][IP2][IP1P1],ggrd[3][IP2P1][IP1P1],ggrd[3][IP2P1][IP1]]
    # print IP1,IP2,IP1P1,IP2P1
    return y,y1,y2,y12


def CMAPSPL(dx,y,n):
    # create spline system
    y2=[0.0] # set lower boundary condition to be "natural"
    u=[0.0]
    dxinv=1.0/dx
    for i in range(1,n-1):
        pinv=1.0/(y2[i-1]+4.0) # should be 2.0? numerical reciept p109
        y2.append(-1.0*pinv)
        u.append(((6*y[i+1]-12*y[i]+6*y[i-1])*dxinv*dxinv-u[i-1])*pinv)
    y2.append(0.0) # set upper boundary condition to be "natural"
    u.append(0.0)
    for i in range(n-2,-1,-1):
        y2[i]=y2[i]*y2[i+1]+u[i]
    return y2

def CMAPSPI(xmin,dx,ya,y2a,x):
    # evaluate the spline and spline derivative
    inx=int(math.floor((x-xmin)/dx))
    a=(xmin+inx*dx+dx-x)/dx
    b=(x-xmin-inx*dx)/dx
    y=a*ya[inx] + b*ya[inx+1] + ( (a*a*a-a)*y2a[inx]+(b*b*b-b)*y2a[inx+1] )*(dx*dx)/6.0
    y1=(ya[inx+1]-ya[inx])/dx - (3*a*a-1)/6*dx*y2a[inx] + (3*b*b-1)/6*dx*y2a[inx+1]
    return y,y1

def setupwt():
    global wt
    wt=  [[   1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
          [   0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], 
          [  -3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1, 0, 0, 0, 0], 
          [   2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0], 
          [   0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
          [   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], 
          [   0, 0, 0, 0,-3, 0, 0, 3, 0, 0, 0, 0,-2, 0, 0,-1], 
          [   0, 0, 0, 0, 2, 0, 0,-2, 0, 0, 0, 0, 1, 0, 0, 1], 
          [  -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
          [   0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0], 
          [   9,-9, 9,-9, 6, 3,-3,-6, 6,-6,-3, 3, 4, 2, 1, 2], 
          [  -6, 6,-6, 6,-4,-2, 2, 4,-3, 3, 3,-3,-2,-1,-1,-2], 
          [   2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 
          [   0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0], 
          [  -6, 6,-6, 6,-3,-3, 3, 3,-4, 4, 2,-2,-2,-2,-1,-1], 
          [   4,-4, 4,-4, 2, 2,-2,-2, 2,-2,-2, 2, 1, 1, 1, 1]]  
    return   

def gcscf(TY,TY1,TY2,TY12,GRES1,GRES2):
    TX=[]
    for i in range(4):
        TX.append(TY[i])
    for i in range(4):
        TX.append(TY1[i]*GRES1)
    for i in range(4):
        TX.append(TY2[i]*GRES2)
    for i in range(4):
        TX.append(TY12[i]*GRES1*GRES2)

    inn=0
    TC=[[],[],[],[]]
    for i in range(4):
        #TC.append([])
	for j in range(4):
            xx=0.0
            a=wt[inn]
	    for k in range(16):
		xx = xx + a[k]*TX[k]
	    inn+=1
	    TC[i].append(xx)

    return TC

def weighting(varcmap,m1,e0):
    nlen=len(varcmap[0])
    # setup CMAP grid
    g1=setcmap2(m1)
    invkT=1.0/(0.0019872041*300)
    w=[]
    for i in range(nframe):
        en=0.0
        for j in range(nlen):
            en = en + ecmapfast(g1,varcmap[i][j][0],varcmap[i][j][1],varcmap[i][j][2],
			    varcmap[i][j][3],varcmap[i][j][4],varcmap[i][j][5])
        w.append(math.exp((e0[i]-en)*invkT))
        #print en
    q=weightingcheck(w)
    return w, q

def energy(varcmap,m1):
    nlen=len(varcmap[0])
    # setup CMAP grid
    g1=setcmap2(m1)
    e=[]
    for i in range(nframe):
        e1=0.0
        for j in range(nlen):
            e1 = e1 + ecmapfast(g1,varcmap[i][j][0],varcmap[i][j][1],varcmap[i][j][2],
			    varcmap[i][j][3],varcmap[i][j][4],varcmap[i][j][5])
        e.append(e1)
    return e

def weightingcheck(a):
# top x% should have less than y% weight
    toppercent=0.05 #0.1,0.01
    cutoff=0.3 #0.4,0.2
    top=int(-1*toppercent*nframe)
    b=sorted(a)[top:]
    c=sum(b)/sum(a)
    #print top,c
    return bool(c<cutoff)

def lefthelical(c1):
# this subroutine readin the phi,psi angles and return an 0/1 array of same size indicating whether
# the residues are left helical or not in each frame
# the definition of left helical: three consecutive residues falling into alpha_L region (30 < phi< 100 and 7 < psi < 67)
    o1=[]
    nlen=len(c1[0])/2 # number of phi,psi data # have to be at least 5
    for c in c1:
        t=[]
        for i in range(nlen):
            if c[2*i]>30 and c[2*i]<100 and c[2*i+1]>7 and c[2*i+1]<67:
                t.append(1)
            else:
                t.append(0)
        tt=[t[0]*t[1]*t[2],t[0]*t[1]*t[2]+t[1]*t[2]*t[3]]
        for i in range(2,nlen-2):
            tt.append(t[i-2]*t[i-1]*t[i]+t[i-1]*t[i]*t[i+1]+t[i]*t[i+1]*t[i+2])
        tt.append(t[nlen-4]*t[nlen-3]*t[nlen-2]+t[nlen-3]*t[nlen-2]*t[nlen-1])
        tt.append(t[nlen-3]*t[nlen-2]*t[nlen-1])
        # we only care about whether there has left helix or not
        ttt=sum(tt)
        if ttt==0:
            o1.append(0)
        else:
            o1.append(1)
    return o1

def lefthelical2(c):
# this subroutine readin the phi,psi angles and return an 0/1 array of same size indicating whether
# the residues are left helical or not in each frame
# the definition of left helical: three consecutive residues falling into alpha_L region (30 < phi< 100 and 7 < psi < 67)
    nlen=len(c)/2 # number of phi,psi data # have to be at least 5
    t=[]
    for i in range(nlen):
        if c[2*i]>30 and c[2*i]<100 and c[2*i+1]>7 and c[2*i+1]<67:
            t.append(1)
        else:
            t.append(0)
    tt=[t[0]*t[1]*t[2],t[0]*t[1]*t[2]+t[1]*t[2]*t[3]]
    for i in range(2,nlen-2):
        tt.append(t[i-2]*t[i-1]*t[i]+t[i-1]*t[i]*t[i+1]+t[i]*t[i+1]*t[i+2])
    tt.append(t[nlen-4]*t[nlen-3]*t[nlen-2]+t[nlen-3]*t[nlen-2]*t[nlen-1])
    tt.append(t[nlen-3]*t[nlen-2]*t[nlen-1])
    # we only care about whether there has left helix or not
    ttt=sum(tt)
    if ttt==0:
        o1=0
    else:
        o1=1
    return o1

def helixdiff(o1,weight):
    t=[]
    tweight=0.0
    for i in range(15):
        t.append(0)

    for ii in range(nframe):
        tweight += weight[ii]
        for i in range(15):
            t[i] += o1[ii][i]*weight[ii]

    for i in range(15):
        t[i] = t[i]/tweight
    tave=sum(t)/15.0
    rms=0.0
    # 300 K experimental data
    texp=[0.194,0.174,0.161,0.262,0.214,0.255,0.245,0.240,0.289,0.230,0.224,0.173,0.141,0.175,0.128,0.000] #Shalongo's 1994 paper
    for i in range(15):
        rms += (t[i]-texp[i])*(t[i]-texp[i])
    a=math.sqrt(rms/15.0)*100.0
    return a,tave

def leftamount(o1,weight):
    tweight=0.0
    oweight=0.0
    for i in range(nframe):
        tweight += weight[i]
        oweight += o1[i]*weight[i]
    return oweight/tweight
        



def rmsd(m1,m2):
# take into two maps (2d lists) and return the square difference between two
    n1=len(m1)
    n2=len(m1[0])
    a=0.0
    for i in range(n1):
        for j in range(n2):
            a += (m1[i][j]-m2[i][j]) * (m1[i][j]-m2[i][j])
    return math.sqrt(a/(n1*n2))

def writeresult(str2,m1):
    for i in range(numspace):
        a=-180+360*i/numspace
        str2.write('!'+str(a)+'\n')
        for j in range(4):
            wline='%8.2f %8.2f %8.2f %8.2f %8.2f\n' % (m1[5*j][i],m1[5*j+1][i],m1[5*j+2][i],m1[5*j+3][i],m1[5*j+4][i])
            str2.write(wline)
        wline='%8.2f %8.2f %8.2f %8.2f\n\n' % (m1[20][i],m1[21][i],m1[22][i],m1[23][i])
        str2.write(wline)

def mcmove1(m1):
# add small random number to each CMAP grid
    n1=len(m1)
    n2=len(m1[0])
    for i in range(n1):
        for j in range(n2):
            m1[i][j] += random.uniform(-0.01, 0.01)
    return

def mcmove2(m1):
# MC move: add a small gaussian to a random CMAP grid 
    n1=len(m1)
    n2=len(m1[0])

    a=100 # preferred region
    b=10 # disfavored region
    c=1  # disallowed region
    # based on Lovell et al: http://kinemage.biochem.duke.edu/downloads/kinfiles/rama/Rama500noGPc.kin
    #  phi -180-165-150-135-120-105 -90 -75 -60 -45 -30 -15  0   15  30  45  60  75  90 105 120 135 150 165
    phipsi=[[a,a,a,a,a,a,a,a,b,b,c,c,c,c,c,c,c,b,c,c,c,c,c,b], # -180
            [b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,b,c,c,c,c,c,c,c], # -165
            [b,b,b,b,b,b,b,b,c,c,c,c,c,c,c,b,b,b,c,c,c,c,c,c], # -150
            [c,c,c,b,b,b,c,c,c,c,c,c,c,c,c,b,b,b,c,c,c,c,c,c], # -135
            [c,c,c,b,b,b,c,c,c,c,c,c,c,c,c,b,b,b,c,c,c,c,c,c], # -120
            [c,c,c,b,b,b,b,c,c,c,c,c,c,c,c,b,b,c,c,c,c,c,c,c], # -105
            [c,c,b,b,b,b,b,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c,c], #  -90
            [c,c,b,b,b,b,b,b,b,b,b,c,c,c,c,c,c,b,c,c,c,c,c,c], #  -75
            [c,c,b,b,b,a,a,a,a,a,b,c,c,c,c,c,b,b,c,c,c,c,c,c], #  -60
            [c,c,b,b,a,a,a,a,a,a,b,c,c,c,c,c,b,b,c,c,c,c,c,c], #  -45
            [c,c,b,b,a,a,a,a,a,a,b,c,c,c,c,c,b,b,c,c,c,c,c,c], #  -30
            [c,c,b,a,a,a,a,a,a,b,b,c,c,c,c,c,b,b,b,c,c,c,c,c], #  -15
            [c,b,b,a,a,a,a,a,a,b,c,c,c,c,c,b,b,b,b,c,c,c,c,c], #    0
            [c,b,b,a,a,a,a,b,b,c,c,c,c,c,c,b,a,a,b,c,c,c,c,c], #   15
            [c,b,b,a,a,a,a,b,c,c,c,c,c,c,b,a,a,a,b,c,c,c,c,c], #   30
            [c,b,b,a,a,b,b,b,b,c,c,c,c,c,b,a,a,b,c,c,c,c,c,c], #   45
            [c,b,b,a,b,b,a,b,b,c,c,c,c,c,b,a,a,b,c,c,c,c,c,c], #   60
            [c,b,a,a,a,b,a,a,b,c,c,c,c,c,b,b,b,b,c,c,c,c,c,c], #   75
            [b,b,a,a,a,a,a,a,b,c,c,c,c,c,c,c,b,c,c,c,c,c,c,c], #   90
            [b,b,a,a,a,a,a,a,b,b,b,c,c,c,c,c,c,c,c,c,c,c,c,c], #  105
            [b,b,a,a,a,a,a,a,a,a,b,c,c,c,c,c,c,c,c,c,c,c,c,c], #  120
            [b,a,a,a,a,a,a,a,a,a,b,c,c,c,c,c,c,c,c,c,c,c,c,c], #  135
            [a,a,a,a,a,a,a,a,a,b,b,c,c,c,c,c,c,b,c,c,c,c,c,b], #  150            
            [a,a,a,a,a,a,a,a,a,b,c,c,c,c,c,c,c,b,c,c,c,c,c,b]] #  165            
            
            
    gheight=random.uniform(-0.5,0.5)
    gweight=random.randint(1,3)

    ii=random.randint(0,n1-1)
    jj=random.randint(0,n2-1)
    #print ii,jj,gheight,gweight
    addgaussian(m1,ii,jj,gheight,gweight)

def mcmove3(m1):
# MC move: add a small gaussian to a trianlge corresponding to the left helix: 30<=phi<=90 && 0<=psi<=60 && phi+psi>=90
    ii=random.randint(12,16) #psi
    jj=random.randint(14,18) #phi
    if ii+jj < 30:
        ii=28-ii
        jj=32-jj

    gheight=random.uniform(-0.2,0.2)
    gweight=random.randint(1,3)

    #print ii,jj,gheight,gweight
    addgaussian(m1,ii,jj,gheight,gweight)

def mcmove4(m1):
# add small fluctuation from point to point in an extended left helix region
    amount=0.01
    namount=-1*amount
    for jj in range(16,19):
        m1[10][jj] += random.uniform(namount, amount)
    for jj in range(16,20):
        m1[11][jj] += random.uniform(namount, amount)
    for jj in range(16,20):
        m1[12][jj] += random.uniform(namount, amount)
    for jj in range(15,21):
        m1[13][jj] += random.uniform(namount, amount)
    for jj in range(14,21):
        m1[14][jj] += random.uniform(namount, amount)
    for jj in range(14,20):
        m1[15][jj] += random.uniform(namount, amount)
    for jj in range(13,20):
        m1[16][jj] += random.uniform(namount, amount)
    for jj in range(12,21):
        m1[17][jj] += random.uniform(namount, amount)
    for jj in range(12,21):
        m1[18][jj] += random.uniform(namount, amount)

def addgaussian(m1,i,j,gh,gw):
# add a gaussian (gh, gw) at [i,j] in a discrete way (assuming f(3sigma)=0)
    n1=len(m1)
    n2=len(m1[0])
    if gw==1:
    # width of 15 degree
        m1[i][j] += gh
    if gw==2:
    # width of 30 degree
        gh10 = gh * 0.325
        gh11 = gh * 0.105
        
        m1[i][j] += gh
        for (a,b) in [(1,0),(-1,0),(0,1),(0,-1)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh10
        for (a,b) in [(1,1),(1,-1),(-1,1),(-1,-1)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh11        
    if gw==3:
    # width of 45 degree
        gh10 = gh * 0.607
        gh11 = gh * 0.368
        gh20 = gh * 0.135
        gh21 = gh * 0.082
        gh22 = gh * 0.018
        
        m1[i][j] += gh
        for (a,b) in [(1,0),(-1,0),(0,1),(0,-1)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh10
        for (a,b) in [(1,1),(1,-1),(-1,1),(-1,-1)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh11        
        for (a,b) in [(2,0),(-2,0),(0,2),(0,-2)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh20
        for (a,b) in [(2,1),(2,-1),(-2,1),(-2,-1),(1,2),(1,-2),(-1,2),(-1,-2)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh21
        for (a,b) in [(2,2),(2,-2),(-2,2),(-2,-2)]:
            ia = (i+a) % n1
            jb = (j+b) % n2
            m1[ia][jb] += gh22

    return        
        

if __name__ == "__main__":

    nresfull=14 # the terminal residues are excluded
    skipres=[2] # Gly and Pro residues with which no standard CMAP is applied
    nres=nresfull-len(skipres)

    global numspace
    numspace=24 # must be even
    wrmsd=2
    setupwt()

    # first step: handling input files: initial CMAP files, coordinate phi,psi files
    try:
        parafile = open(sys.argv[1],'r')
    except IOError:                     #catch exception
        print ('Files do not exist!\n')
    try:
        coorfile = open(sys.argv[2],'r')
    except IOError:                     #catch exception
        print ('Files do not exist!\n')
    try:
        writefile = open(sys.argv[3],'w')
    except IOError:                     #catch exception
        print ('Files do not exist!\n')

    # read coordinates files (with respect to energy evaluation)
    global nframe
    coor1=[]
    coor2=[]
    ii=0
    coor=[]
    obs=[]
    for line in coorfile.xreadlines():
        coor3=map(float,string.split(line))
        coor1.append(coor3[0])
        coor1.append(coor3[1])
        if ii not in skipres:
	    coor2.append(convertvar(coor3[0],coor3[1]))
        # update the index
        if ii < nresfull-1:
            ii += 1
        elif ii == nresfull-1:
            obs.append(lefthelical2(coor1))
            coor.append(coor2)
            coor1=[]
            coor2=[]
            ii=0
        else:
            print ('Error in reading coordinates!\n')        
    nframe=len(obs)        

    # set up initial CMAP 
    map_0=readcmap(parafile)


    # initialize
    estart=energy(coor,map_0)

    map_old=copy.deepcopy(map_0)  # this is what we suppose to do
    t1 = sum(obs)/float(nframe)
    tar_old = math.log(t1)
    tar_new = tar_old
    tar_best = tar_old
    prob=0
    accepted=1

    ###### do the monte carlo simulated annealing
    tempr0 = 10.0
    nstep = 50000
    #print some basic information
    print sys.argv[2], nframe, wrmsd
    print tempr0, nstep

    print "%s" % ( "Starting Monte Carlo fitting" )
    print "%11s%11s%11s%11s%11s%11s" % ( "MC step", "tempr", "p", "accepted?", "RMSE", "RMSE_best" )
    step = 0
    while step < nstep:
        tempr = tempr0 * math.exp(-1.0* ( float( step ) / ( float( nstep ) / 4.0 ) ))
        map_new=copy.deepcopy(map_old)
        mcmove4(map_new)
        weight_new,qweight=weighting(coor,map_new,estart)
        if qweight == True:
            t1 = leftamount(obs,weight_new)
            t3 = rmsd(map_new,map_0)
            tar_new=math.log(t1)+wrmsd*t3
            #print tempr, t1, t2, t3, qweight        
            dtar= tar_new - tar_old
            if step == 0:
                prob  = 0.0  # required to prevent overflow
            else:
                dtar= tar_new - tar_old
                boltz = -1.0 * dtar / ( 0.001987 * tempr )
                prob  = math.exp( boltz )
            prob0 = random.uniform(0.0,1.0)
            #print "%11i%10.5f%10.5f%10.5f%10.5f%10.5f" % (step, prob, prob0, tar_new, tar_old, dtar)
            if dtar < 0.0:
                accepted = 1
                prob = 1.0
            elif prob0 < prob:
                accepted = 1
            else:
                accepted = 0
            # print accepted
            if accepted:
                tar_old = tar_new
                map_old=copy.deepcopy(map_new) # memory inefficient?
                if tar_new < tar_best:
                    tar_best = tar_new
                    map_best=copy.deepcopy(map_new)
        if ( step % 10 ) == 0:
            print "%11i%11.1f%11.4f%11i%11.2f%11.2f%11.4f%11.4f" % ( step, tempr, prob, accepted, tar_new, tar_best, t1, t3)
        step = step + 1
    print "%11i%11.1f%11.4f%11i%11.2f%11.2f%11.4f%11.4f" % ( step, tempr, prob, accepted, tar_new, tar_best, t1, t3)


    # now writeout the optimized CMAP
    writeresult(writefile,map_best)
    writefile.write('\nFurther Info:\n')
    wline='based on '+str(nframe)+' frames and w='+str(wrmsd)
    writefile.write(wline)
    wline='\nRMSD with Init CMAP: '+str(rmsd(map_best,map_0))
    writefile.write(wline)
    weight_new,qweight=weighting(coor,map_best,estart)
    t1=leftamount(obs,weight_new)
    wline='\nPrediction of left helical contents: '+str(t1)
    writefile.write(wline)

