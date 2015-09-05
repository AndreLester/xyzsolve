#!/usr/bin/env python
#Adapted from fortran dorpri5 package

import math
#Dormand-Prince Runge-Kutta integrator
from Limits import *

CDOPRI_C2=0.20
CDOPRI_C3=0.30
CDOPRI_C4=0.80
CDOPRI_C5=8.0/9.0
CDOPRI_A21=0.20
CDOPRI_A31=3.0/40.0
CDOPRI_A32=9.0/40.0
CDOPRI_A41=44.0/45.0
CDOPRI_A42=-56.0/15.0
CDOPRI_A43=32.0/9.0
CDOPRI_A51=19372.0/6561.0
CDOPRI_A52=-25360.0/2187.0
CDOPRI_A53=64448.0/6561.0
CDOPRI_A54=-212.0/729.0
CDOPRI_A61=9017.0/3168.0
CDOPRI_A62=-355.0/33.0
CDOPRI_A63=46732.0/5247.0
CDOPRI_A64=49.0/176.0
CDOPRI_A65=-5103.0/18656.0
CDOPRI_A71=35.0/384.0
CDOPRI_A73=500.0/1113.0
CDOPRI_A74=125.0/192.0
CDOPRI_A75=-2187.0/6784.0
CDOPRI_A76=11.0/84.0
CDOPRI_E1=71.0/57600.0
CDOPRI_E3=-71.0/16695.0
CDOPRI_E4=71.0/1920.0
CDOPRI_E5=-17253.0/339200.0
CDOPRI_E6=22.0/525.0
CDOPRI_E7=-1.0/40.0
CDOPRI_D1=-12715105075.0/11282082432.0
CDOPRI_D3=87487479700.0/32700410799.0
CDOPRI_D4=-10690763975.0/1880347072.0
CDOPRI_D5=701980252875.0/199316789632.0
CDOPRI_D6=-1453857185.0/822651844.0
CDOPRI_D7=69997945.0/29380423.0



def sign(x, y):
    if y >= 0.0:
        return abs(x)
    else:
        return -abs(x)



"""
Default Parameters:
NMAX = 100000 (Maximum number of function calls)
SAFE = 0.9 (Bounds 1.0e-4 to 1.0)
FAC1 = 0.2 (No bounds available)
FAC2 = 10.0 (No bounds available)
BETA = 0.04 (Bounds <0.2)
HMAX = XEND-X (Maximum step size)
NSTIFF = 1000 (Iterations before equation set is considered stiff)
"""

"""
Notes:
When ITOL is 0, RTOL and ATOL are numbers
When ITOL is not 0, RTOL and ATOL are vectors containing N numbers

This info is really not useful to the user
------------------------------------------
NFCN is number of function calls
NSTEP is number of steps
NACCPT is number of accepted steps
NREJCT is number of rejected steps

"""

def DOPCOR(N, FCN, X, Y, XEND, HMAX, H, RTOL, ATOL, ITOL, IPRINT,
           SOLOUT, NMAX, NSTIFF, SAFE, BETA, FAC1, FAC2, 
           CONT, ICOMP, NRD, RPAR, NFCN, NSTEP, NACCPT, NREJCT):
    
    UROUND = EPSMCH
    FACOLD = 1.E-4  
    EXPO1 = 0.2 - BETA*0.75
    FACC1 = 1.0/FAC1
    FACC2 = 1.0/FAC2
    POSNEG = sign(1.0, XEND-X)
    IDID = 0
    
    IOUT = 0
    if SOLOUT is not None:
        IOUT = 1
    
    Y1 = [0.0]*N
    K1 = [0.0]*N
    K2 = [0.0]*N
    K3 = [0.0]*N
    K4 = [0.0]*N
    K5 = [0.0]*N
    K6 = [0.0]*N
    YSTI = [0.0]*N
    
    if ITOL > 0:
        ATOLI = ATOL[0]
        RTOLI = RTOL[0]    
    else:
        ATOLI = ATOL
        RTOLI = RTOL
    
    LAST = False 
    HLAMB = 0.0
    IASTI = 0
    FCN(X, Y, K1, RPAR)
    HMAX = abs(HMAX)     
    IORD = 5
    
    if (H == 0.0):
        H = HINIT(N, FCN, X, Y, XEND, POSNEG, K1, K2, K3, IORD, HMAX,
                  ATOL, RTOL, ITOL, RPAR)
    
    NFCN = NFCN+2
    REJECT = False
    XOLD = X
    
    if (IOUT != 0): 
        IRTRN = 1
        HOUT = H
        SOLOUT(NACCPT+1, XOLD, X, Y, N, CONT, ICOMP, NRD, RPAR, IPAR, IRTRN)
        if (IRTRN < 0):
            IDID = 2
            return IDID
    else:
        IRTRN = 0 
    
    LOOP = True
    
    while(LOOP):
        IRTRN = 2
        if (NSTEP > NMAX):
            print 'MORE THAN NMAX =',NMAX,' STEPS ARE NEEDED'
            IDID = -2
            return IDID
        
        if (0.10*abs(H) <= abs(X)*UROUND):
            print 'STEP SIZE TOO SMALL, H=', H
            IDID = -3
            return IDID
        
        if ((X + 1.01*H - XEND)*POSNEG > 0.0):
            H = XEND-X 
            LAST = True
        
        NSTEP = NSTEP + 1
        
# --- THE FIRST 6 STAGES
        if (IRTRN >= 2):
            FCN(X, Y, K1, RPAR)
        
        for I in xrange(0, N):
            Y1[I] = Y[I]+ H*CDOPRI_A21*K1[I]
        FCN(X+CDOPRI_C2*H, Y1, K2, RPAR)
        
        for I in xrange(0, N):
            Y1[I] = Y[I]+H*(CDOPRI_A31*K1[I]+CDOPRI_A32*K2[I])
        FCN(X+CDOPRI_C3*H,Y1,K3,RPAR)
        
        for I in xrange(0, N):
            Y1[I] = Y[I]+H*(CDOPRI_A41*K1[I]+CDOPRI_A42*K2[I]+CDOPRI_A43*K3[I])        
        FCN(X+CDOPRI_C4*H, Y1, K4, RPAR)
        
        for I in xrange(0, N):
            Y1[I] = Y[I]+H*(CDOPRI_A51*K1[I]+CDOPRI_A52*K2[I]+CDOPRI_A53*K3[I]+CDOPRI_A54*K4[I])
        FCN(X+CDOPRI_C5*H, Y1, K5, RPAR) 
        
        for I in xrange(0, N):
            YSTI[I] = Y[I]+H*(CDOPRI_A61*K1[I]+CDOPRI_A62*K2[I]+CDOPRI_A63*K3[I]+CDOPRI_A64*K4[I]+CDOPRI_A65*K5[I])
        XPH = X+H
        FCN(XPH, YSTI, K6, RPAR)
        
        for I in xrange(0, N):
            Y1[I] = Y[I]+H*(CDOPRI_A71*K1[I]+CDOPRI_A73*K3[I]+CDOPRI_A74*K4[I]+CDOPRI_A75*K5[I]+CDOPRI_A76*K6[I])
        FCN(XPH, Y1, K2, RPAR)
        
        if (IOUT >= 2):
            for J in xrange(1, NRD+1):
                I = ICOMP[J-1]
                CONT[(4*NRD+J)-1]=H*(CDOPRI_D1*K1[I-1]+CDOPRI_D3*K3[I-1]+CDOPRI_D4*K4[I-1]+CDOPRI_D5*K5[I-1]+CDOPRI_D6*K6[I-1]+CDOPRI_D7*K2[I-1])
                
        
        for I in xrange(0, N):
            K4[I]=(CDOPRI_E1*K1[I]+CDOPRI_E3*K3[I]+CDOPRI_E4*K4[I]+CDOPRI_E5*K5[I]+CDOPRI_E6*K6[I]+CDOPRI_E7*K2[I])*H
        
        NFCN = NFCN + 6
        
        ERR = 0.0
        if (ITOL == 0.0):
            for I in xrange(0, N):
                SK = ATOLI + RTOLI*max(abs(Y[I]), abs(Y1[I]))
                ERR += (K4[I]/SK)*(K4[I]/SK)
        else:
            for I in xrange(0, N):
                SK = ATOL[I] + RTOL[I]*max(abs(Y[I]), abs(Y1[I]))
                ERR += (K4[I]/SK)*(K4[I]/SK)
        
        ERR = math.sqrt(ERR/N)
        
        FAC11 = math.pow(ERR, EXPO1)
        
        FAC = FAC11/(math.pow(FACOLD, BETA))
        
        FAC = max(FACC2, min(FACC1,FAC/SAFE))
        HNEW = H/FAC
        
        if (ERR <= 1.0):
            
            FACOLD = max(ERR, 1.0e-4)
            NACCPT=NACCPT + 1
            
            if (NACCPT % NSTIFF == 0 or IASTI>0):
                STNUM = 0.0
                STDEN = 0.0
            
                for I in xrange(0, N): 
                    STNUM = STNUM+(K2[I]-K6[I])*(K2[I]-K6[I])
                    STDEN = STDEN+(Y1[I]-YSTI[I])*(Y1[I]-YSTI[I])
                
                if (STDEN > 0.00):
                    HLAMB = H*math.sqrt(STNUM/STDEN) 
                
                if (HLAMB > 3.25):
                    NONSTI = 0
                    IASTI = IASTI + 1  
                    if (IASTI == 15):
                        if (IPRINT > 0):
                            print 'THE PROBLEM SEEMS TO BECOME STIFF AT X = ', X
                        if (IPRINT <= 0):
                            IDID = -4
                            return IDID
                else:
                    NONSTI = NONSTI + 1  
                    if (NONSTI == 6):
                        IASTI=0
                
            if (IOUT >= 2):
                for J in xrange(1, NRD+1):
                    I = ICOMP[J-1]
                    YD0 = Y[I]
                    YDIFF = Y1[I] - YD0
                    BSPL = H*K1[I] - YDIFF 
                    CONT[J-1] = Y[I]
                    CONT[NRD+J-1] = YDIFF
                    CONT[2*NRD+J-1] = BSPL
                    CONT[3*NRD+J-1] = -H*K2(I)+YDIFF-BSPL
            
            for I in xrange(0, N):
                K1[I] = K2[I]
                Y[I] = Y1[I]
            
            XOLD = X
            X = XPH
            
            if (IOUT != 0):
                HOUT = H
                SOLOUT(NACCPT+1,XOLD,X,Y,N,CONT,ICOMP,NRD, RPAR,IPAR,IRTRN)
                if (IRTRN < 0):
                    IDID = 2
                    return IDID
            
            if (LAST):
                H = HNEW
                IDID = 1
                return
            
            if (abs(HNEW) > HMAX):
                HNEW = POSNEG*HMAX
            
            if (REJECT):
                HNEW = POSNEG*min(abs(HNEW), abs(H))
            
            REJECT = False
            
        else:
            HNEW = H/min(FACC1,FAC11/SAFE)
            REJECT = True
            if (NACCPT >= 1):
                NREJCT = NREJCT + 1
            LAST = False
        
        H = HNEW
    
    #Return to the main loop
    
    IDID = -2
    return IDID
 


def HINIT(N, FCN, X, Y, XEND, POSNEG, F0, F1, Y1, IORD, HMAX, 
          ATOL, RTOL, ITOL, RPAR):
        
    DNF=0.0
    DNY=0.0
    if ITOL > 0:
        ATOLI = ATOL[0]
        RTOLI = RTOL[0]
    else:
        ATOLI = ATOL
        RTOLI = RTOL

    if (ITOL == 0):   
        for I in xrange(0, N): 
            SK = ATOLI + RTOLI*abs(Y[I])
            DNF = DNF + (F0[I]/SK)*(F0[I]/SK)
            DNY = DNY + (Y[I]/SK)*(Y[I]/SK)
    else:
        for I in xrange(0, N): 
            SK = ATOL[I]+RTOL[I]*abs(Y[I])
            DNF = DNF+(F0[I]/SK)*(F0[I]/SK)
            DNY = DNY+(Y[I]/SK)*(Y[I]/SK)
            
    if (DNF <= 1.0E-10 or DNY <= 1.0E-10):
        H = 1.0E-6
    else:
        H = math.sqrt(DNY/DNF)*0.01 
    
    H = min(H, HMAX)
    H = sign(H, POSNEG) 
    
    for I in xrange(0, N):
        Y1[I] = Y[I] + H*F0[I]
    
    FCN(X+H, Y1, F1, RPAR)

    DER2 = 0.0
    if (ITOL == 0):   
        for I in xrange(0, N): 
            SK = ATOLI + RTOLI*abs(Y[I])
            DER2 = DER2 + ((F1[I]-F0[I])/SK)*((F1[I]-F0[I])/SK)   
    else:
        for I in xrange(0, N): 
            SK = ATOL[I] + RTOL[I]*abs(Y[I])
            DER2 = DER2 + ((F1[I]-F0[I])/SK)*((F1[I]-F0[I])/SK)
    
    DER2 = math.sqrt(DER2)/H
    
    DER12 = max(abs(DER2), math.sqrt(DNF))
    if (DER12 <= 1.0e-15):
        H1 = max(1.0e-6, abs(H)*1.0e-3)
    else:
        H1=math.pow((0.01/DER12),(1.00/IORD))
        
    H = min(100.0*abs(H), H1, HMAX)
    
    return sign(H, POSNEG)


def CONTD5(II, X, CON, ICOMP, ND, XOLD, H):
    I = 0
    for J in xrange(0, ND):
        if ICOMP[J] == II:
            I = J
    if (I == 0):
        print 'NO DENSE OUTPUT AVAILABLE FOR COMP. ', II
        return
    
    THETA=(X-XOLD)/H
    THETA1 = 1.0 - THETA
    return CON[I] + THETA*(CON[ND+I-1]+THETA1*(CON[2*ND+I-1]+THETA*(CON[3*ND+I-1]+THETA1*CON[4*ND+I-1])))


def dormand_prince_odeint(func, tstart, tstop, y0, step=0.0, par=None):
    """
    Default Parameters:
    NMAX = 100000 (Maximum number of function calls)
    SAFE = 0.9 (Bounds 1.0e-4 to 1.0)
    FAC1 = 0.2 (No bounds available)
    FAC2 = 10.0 (No bounds available)
    BETA = 0.04 (Bounds <0.2)
    HMAX = XEND-X (Maximum step size)
    NSTIFF = 1000 (Iterations before equation set is considered stiff)
    """

    """
    Notes:
    When ITOL is 0, RTOL and ATOL are numbers
    When ITOL is not 0, RTOL and ATOL are vectors containing N numbers

    This info is really not useful to the user
    ------------------------------------------
    NFCN is number of function calls
    NSTEP is number of steps
    NACCPT is number of accepted steps
    NREJCT is number of rejected steps

    """
    n = len(y0)
    NMAX = 100000
    SAFE = 0.9
    FAC1 = 0.2
    FAC2 = 10.0
    BETA = 0.04
    H = step
    HMAX = step
    if HMAX == 0.0:
        HMAX = tstop - tstart
    NSTIFF = 1000
    NFCN = 0
    NSTEP = 0
    NACCPT = 0
    NREJCT = 0
    ITOL = 0
    IPRINT = 1
    RTOL = 1e-4
    ATOL = 1e-6
    NRD = 0
    ICOMP = 0
    CONT = 0
    return DOPCOR(n, func, tstart, y0, tstop, HMAX, H, RTOL, ATOL, ITOL, IPRINT,
                  None, NMAX, NSTIFF, SAFE, BETA, FAC1, FAC2, 
                  CONT, ICOMP, NRD, par, NFCN, NSTEP, NACCPT, NREJCT)

def step(t, limit):
    result = 0.0
    if t == limit:
        result = 0.5
    if t > limit:
        result = 1.0 

    return result

def spike(t, limit):
    if t == limit:
        return 1.0
    
    return 0.0

def test():
    y0 = [0.994, 0.0, 0.0, -2.00158510637908252240537862224]

    def func(t, y, f, par):
        amu = 0.012277471
        amup = 1.0 - amu
        
        r1 = ((y[0] + amu)**2 + y[1]**2 ) ** (3.0/2.0)
        r2 = ((y[0] - amup)**2 + y[1]**2 ) ** (3.0/2.0)
        
        
        f[0] = y[2]
        f[1] = y[3]
        f[2] = y[0] + 2.0 * y[3] - amup * (y[0]+amu) / r1 - amu * (y[0]-amup) / r2
        f[3] = y[1] - 2.0 * y[2] - amup * y[1] / r1 - amu * y[1] / r2
        
        print t,'\t'.join([str(y_value) for y_value in y])
        

    dormand_prince_odeint(func, 0, 17.0652165601579625588917206249, y0)


def test2():
    y = [1.0]
    def func(t, y, f, par):
        #print t,'\t'.join([str(y_value) for y_value in y])
        f[0] = -y[0] + math.exp(1/y[0]) - step(t, 5) + (step(t, 10))*y[0]
    
    delta_t = 0.5
    t = 0.0
    max_t = 20.0
    while t <= max_t:
        dormand_prince_odeint(func, t, t+delta_t, y)
        print t,'\t'.join([str(y_value) for y_value in y])
        t += delta_t

    
if __name__ == '__main__':
    test2()

