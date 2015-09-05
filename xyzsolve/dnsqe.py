#!/usr/bin/env python
# Adapted from the dnsqe fortran package

import math
from Limits import *



def frange(start, stop, inc=1):
    if start <= stop:
        return xrange(start, stop+1, inc)
    else:
        if inc == 1:
            inc = -inc;
        return xrange(start, stop-1, inc)

def dqform (m, n,  q, ldq, wa):    
    minmn = min(m, n)
    if (minmn >=2):
        for j in xrange(2, minmn+1):
            jm1 = j - 1
            for i in xrange(1, jm1+1):
                q[i-1][j-1] = 0.0
        

    np1 = n + 1
    if (m >= np1):
        for j in xrange(np1,m+1):
            for i in xrange(1, m+1):
                q[i-1][j-1] = 0.0
            
            q[j-1][j-1] = 0.0
        

    for l in xrange(1, minmn+1):
        k = minmn - l + 1
        for i in xrange(k, m+1):
            wa[i-1] = q[i-1][k-1]
            q[i-1][k-1] = 0.0
         
        q[k-1][k-1] = 1.0
        if (wa[k-1] != 0.0): 
            for j in xrange(k, m+1):
                sum_dqform = 0.0
                for i in xrange(k, m+1):
                    sum_dqform += q[i-1][j-1]*wa[i-1]
                 
                temp = sum_dqform/wa[k-1]
                for i in xrange(k, m+1):
                    q[i-1][j-1] = q[i-1][j-1] - temp*wa[i-1]
                 
             
    return

def column_denorm(row_count, a, start_row, column):
    sum_column_denorm = 0.0
    for i in xrange(1, row_count+1):
        sum_column_denorm += a[(start_row+i-1)-1][column-1]*a[(start_row+i-1)-1][column-1]
    
    return math.sqrt(sum_column_denorm)


def vector_denorm(n, a):
    sum_vector_denorm = 0.0
    for i in xrange(1, n+1):
        sum_vector_denorm += a[i-1]*a[i-1]
    
    return math.sqrt(sum_vector_denorm)

def dqrfrac (m, n, a, lda, pivot, ipvt, lipvt, sigma, acnorm, wa):
            
    minmn = 0
    kmax = 0.0
    
    for j in xrange(1, n+1):
        acnorm[j-1] = column_denorm(m, a, 1, j)
        sigma[j-1] = acnorm[j-1]
        wa[j-1] = sigma[j-1]
        if (pivot):
            ipvt[j-1] = j
        
    
    minmn = min(m, n)
    for j in xrange(1, minmn+1):
        if (pivot):
            kmax = j
            for k in xrange(j, n+1):
                if (sigma[k-1] > sigma[kmax-1]):
                    kmax = k;
                 
             
            if (kmax != j):
                k = 0
                temp = 0.0
                for i in xrange(1, m+1):
                    temp = a[i-1][j-1]
                    a[i-1][j-1] = a[i-1][kmax-1]
                    a[i-1][kmax-1] = temp
                 

                sigma[kmax-1] = sigma[j-1]
                wa[kmax-1] = wa[j-1]
                k = ipvt[j-1]
                ipvt[j-1] = ipvt[kmax-1]
                ipvt[kmax-1] = k
             
         
        ajnorm = column_denorm(m-j+1, a, j, j)
        if (ajnorm != 0.0):
            if (a[j-1][j-1] < 0.0):
                ajnorm = -ajnorm
             
            for i in xrange(j, m+1):
                a[i-1][j-1] = a[i-1][j-1]/ajnorm
             
            a[j-1][j-1] = a[j-1][j-1] + 1.0

            jp1 = j + 1
            if (n >= jp1):
                for k in xrange(jp1, n+1):
                    sum_ = 0.0
                    for i in xrange(j, m+1):
                        sum_ += a[i-1][j-1]*a[i-1][k-1]
                    
                    temp = sum_/a[j-1][j-1]
                    for i in xrange(j, m+1):
                        a[i-1][k-1] = a[i-1][k-1] - temp*a[i-1][j-1]
                     
                    if (pivot and sigma[k-1]!= 0.0):
                        temp = a[j-1][k-1]/sigma[k-1]
                        sigma[k-1] = sigma[k-1]*math.sqrt(max(0.0, 1.0 - temp*temp))
                        if (0.05*(sigma[k-1]/wa[k-1])*(sigma[k-1]/wa[k-1]) <= EPSMCH):
                            sigma[k-1] = column_denorm(m-j, a, jp1, k)
                            wa[k-1] = sigma[k-1]
                         
                     
        sigma[j-1] = -ajnorm
     
    return

def d1updt (m, n, s, ls, u, v, w, sing):
    
    jj = int(math.floor((n*(2*m - n + 1))/2.0) - (m - n))
    l = jj

    for i in xrange(n, m+1):
        w[i-1] = s[l-1]
        l = l + 1
     

    nm1 = n - 1
    if (nm1 >= 1):

        for nmj in xrange(1, nm1+1):
            j = n - nmj
            jj = jj - (m - j + 1)
            w[j-1] = 0

            if (v[j-1] != 0.0):
                cotan = 0.0
                tan = 0.0
                sin = 0.0
                cos = 0.0
                tau = 0.0

                if (abs(v[n-1]) < abs(v[j-1])):
                    cotan = v[n-1]/v[j-1]
                    sin = 0.5/math.sqrt(0.25 + 0.25*(cotan*cotan))
                    cos = sin*cotan
                    tau = 1.0
                    if (abs(cos)*GIANT > 1):
                        tau = 1.0/cos
                     
                 
                else:
                    tan = v[j-1]/v[n-1]
                    cos = 0.5/math.sqrt(0.25 + 0.25*(tan*tan))
                    sin = cos*tan
                    tau = sin
                 

                v[n-1] = sin*v[j-1] + cos*v[n-1]
                v[j-1] = tau

                l = jj
                for i in xrange(j, m+1):
                    temp = cos*s[l-1] - sin*w[i-1]
                    w[i-1] = sin*s[l-1] + cos*w[i-1]
                    s[l-1] = temp
                    l = l + 1
                 
             
         
     
    for i in xrange(1, m+1):
        w[i-1] = w[i-1] + v[n-1]*u[i-1]
    

    sing[0] = False
    if (nm1 >= 1):
        for j in xrange(1, nm1+1):
            if (w[j-1] != 0.0):
                cotan = 0.0
                sin = 0.0
                cos = 0.0
                tan = 0.0 
                tau = 0.0 
                temp = 0.0 
                if (abs(s[jj-1]) < abs(w[j-1])):
                    cotan = s[jj-1]/w[j-1]
                    sin = 0.5/math.sqrt(0.25 + 0.25*(cotan*cotan))
                    cos = sin*cotan
                    if (abs(cos)*GIANT > 1):
                        tau = 1.0/cos
                     
                else:
                    tan = w[j-1]/s[jj-1]
                    cos = 0.5/math.sqrt(0.25 + 0.25*(tan*tan))
                    sin = cos*tan
                    tau = sin
                 

                l = jj
                for i in xrange(j, m+1):
                    temp = cos*s[l-1] + sin*w[i-1]
                    w[i-1] = -sin*s[l-1] + cos*w[i-1]
                    s[l-1] = temp
                    l = l + 1
                 
                w[j-1] = tau 
             
            if (s[jj-1] == 0.0):
                sing[0] = True
            
            jj = jj + (m - j + 1) 
         
     
    l = jj
    for i in xrange(n, m+1):
        s[l-1] = w[i-1]
        l = l + 1
     
    if (s[jj-1] == 0.0):
        sing[0] = True
    
    return


def d1mpyq_m (m, n, a, lda, v, w):
    nm1 = n - 1
    j = 0
    cos = 0.0
    sin = 0.0
    temp = 0.0

    if (nm1 >= 1):
        for nmj in xrange(1, nm1+1):
            j = n - nmj
            if (abs(v[j-1]) > 1.0):
                cos = 1.0/v[j-1]
             
            if (abs(v[j-1]) > 1):
                sin = math.sqrt(1.0 - cos*cos)
             
            if (abs(v[j-1]) < 1):
                sin = v[j-1]
             
            if (abs(v[j-1]) < 1):
                cos = math.sqrt(1.0 - sin*sin);
             

            for i in xrange(1, m+1):
                temp = cos*a[i-1][j-1] - sin*a[i-1][n-1]
                a[i-1][n-1] = sin*a[i-1][j-1] + cos*a[i-1][n-1]
                a[i-1][j-1] = temp
             

        for j in xrange(1, nm1+1):
            if (abs(w[j-1]) > 1.0):
                cos = 1.0/w[j-1]
             
            if (abs(w[j-1]) > 1.0):
                sin = math.sqrt(1.0 - cos*cos)
             
            if (abs(w[j-1]) < 1.0):
                sin = w[j-1]
             
            if (abs(w[j-1]) < 1.0):
                cos = math.sqrt(1.0 - sin*sin)
             
            for i in xrange(1, m+1):
                temp = cos*a[i-1][j-1] + sin*a[i-1][n-1]
                a[i-1][n-1] = -sin*a[i-1][j-1] + cos*a[i-1][n-1]
                a[i-1][j-1] = temp
             
         
    return 
 

def d1mpyq (m, n, a, lda, v, w):
    nm1 = n - 1
    j = 0
    cos = 0.0
    sin = 0.0
    temp = 0.0

    if (nm1 >= 1):
        for nmj in xrange(1, nm1+1):
            j = n - nmj
            if (abs(v[j-1]) > 1.0):
                cos = 1.0/v[j-1]
             
            if (abs(v[j-1]) > 1):
                sin = math.sqrt(1.0 - cos*cos)
             
            if (abs(v[j-1]) < 1):
                sin = v[j-1]
             
            if (abs(v[j-1]) < 1):
                cos = math.sqrt(1.0 - sin*sin)
             

            for i in xrange(1, m+1):
                temp = cos*a[j-1] - sin*a[n-1]
                a[n-1] = sin*a[j-1] + cos*a[n-1]
                a[j-1] = temp
             
         

        for j in xrange(1, nm1+1):
            if (abs(w[j-1]) > 1.0):
                cos = 1.0/w[j-1]
             
            if (abs(w[j-1]) > 1.0):
                sin = math.sqrt(1.0 - cos*cos)
             
            if (abs(w[j-1]) < 1.0):
                sin = w[j-1]
             
            if (abs(w[j-1]) < 1.0):
                cos = math.sqrt(1.0 - sin*sin)
             
            for i in xrange(1, m+1):
                temp = cos*a[j-1] + sin*a[n-1]
                a[n-1] = -sin*a[j-1] + cos*a[n-1]
                a[j-1] = temp
             
    
    return
    


def ddoglg (n, r, lr, diag, qtb, delta, x, wa1, wa2):
    
    jj = int(math.floor((n*(n + 1))/2) + 1)
    j = 0
    jp1 = 0
    l = 0
    sum_ = 0.0
    temp = 0.0
    qnorm = 0.0
    gnorm = 0.0
    sgnorm = 0.0
    bnorm = 0.0
    alpha = 0.0

    for k in xrange(1, n+1):
        j = n - k + 1
        jp1 = j + 1
        jj = jj - k
        l = jj + 1
        sum_ = 0.0
        if (n >= jp1):
            for i in xrange(jp1, n+1):
                sum_ += r[l-1]*x[i-1];
                l += 1;
             
         
        temp = r[jj-1]
        if (temp == 0.0):
            l = j
            for i in xrange(1, j+1):
                temp = max(temp, abs(r[l-1]))
                l = l + n - i
             
            temp = EPSMCH*temp
            if (temp == 0.0):
                temp = EPSMCH
             
         
        x[j-1] = (qtb[j-1] - sum_)/temp
     
    for j in xrange(1, n+1):
        wa1[j-1] = 0.0
        wa2[j-1] = diag[j-1] * x[j-1]
     
    qnorm = vector_denorm(n, wa2)
    if (qnorm > delta):

        l = 1
        for j in xrange(1, n+1):
            temp = qtb[j-1]
            for i in xrange(j, n+1):
                wa1[i-1] = wa1[i-1] + r[l-1]*temp
                l = l + 1;
            
            wa1[j-1] = wa1[j-1]/diag[j-1]
         

        gnorm = vector_denorm(n, wa1)
        sgnorm = 0.0
        alpha = delta/qnorm

        if (gnorm != 0.0):

            for j in xrange(1, n+1):
                wa1[j-1] = (wa1[j-1]/gnorm)/diag[j-1]
             

            l = 1
            for j in xrange(1, n+1):
                sum_ = 0.0
                for i in xrange(j, n+1):
                    sum_ += r[l-1]*wa1[i-1]
                    l = l + 1
                 
                wa2[j-1] = sum_
             

            temp = vector_denorm(n, wa2)
            sgnorm = (gnorm/temp)/temp

            alpha = 0.0
            if (sgnorm < delta):
                bnorm = vector_denorm(n, qtb)
                temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
                temp = temp - (delta/qnorm)*((sgnorm/delta)*(sgnorm/delta))
                temp = temp + math.sqrt( (temp - (delta/qnorm))*(temp - (delta/qnorm)) + (1 - (delta/qnorm)*(delta/qnorm))*(1 - (sgnorm/delta)*(sgnorm/delta)))
                alpha = ((delta/qnorm)*(1 - (sgnorm/delta)*(sgnorm/delta)))
             
         

        temp = (1.0 - alpha)*min(sgnorm, delta)
        for j in xrange(1, n+1):
            x[j-1] = temp*wa1[j-1] + alpha*x[j-1]
        
    return



def dfdjc1 (fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu, epsfcn, wa1, wa2):
    
    eps = math.sqrt(max(epsfcn, EPSMCH))
    msum = ml + mu + 1
    temp = 0.0
    h = 0.0
    
    if (msum >= n):
        for j in xrange(1, n+1):
            temp = x[j-1]
            h = eps*abs(temp)
            if (h == 0.0):
                h = eps
            
            x[j-1] = temp + h
            fcn (n, x, wa1, iflag)
            if (iflag[0] < 0):
                return
            
            x[j-1] = temp
            for i in xrange(1, n+1):
                fjac[i-1][j-1] = (wa1[i-1] -fvec[i-1])/h
             
        
    else:
        for k in xrange(1, msum+1):
            for j in xrange(k, n+1):
                wa2[j-1] = x[j-1]
                h = eps*abs(wa2[j-1])
                if (h == 0.0):
                    h = eps
                 
                x[j-1] = wa2[j-1] + h
             
            fcn (n, x, wa1, iflag) 
            if (iflag[0] < 0):
                return;
              
            for j in xrange(k, n+1):
                x[j-1] = wa2[j-1]
                h = eps*abs(wa2[j-1])
                if (h == 0.0):
                    h = eps
                 
                for i in xrange(1, n+1):
                    fjac[i-1][j-1] = 0.0
                    if (i >= j-mu and i<= j + ml):
                        fjac[i-1][j-1] = (wa1[i-1] - fvec[i-1])/h;
                     
                
    return


def dnsq(fcn, jac, iopt, n, x, fvec, fjac, ldfjac,
         xtol, maxfev, ml, mu, epsfcn, diag, mode, factor, 
         nprint, info, nfev, njev, r, lr, qtf, wa1, wa2, wa3, wa4, bounds):
    
    iflag = [0]
    outer_loop_valid = True
    
    #int i, j, l, jm1;
    #double fnorm = 0.0, fnorm1, sum_, delta = 0.0, xnorm = 0.0, temp;
    #double actred, prered, ratio;
    
    info[0] = 0;
    iflag[0] = 0;
    nfev[0] = 0;
    njev[0] = 0;

    if (iopt < 1  or iopt > 2  or n <= 0):
        outer_loop_valid = False
    elif (xtol < 0 or maxfev <= 0 or ml < 0):
        outer_loop_valid = False
    elif (ml < 0 or factor <= 0 or ldfjac < n):
        outer_loop_valid = False
    elif (lr < int(math.floor((n*(n+1))/2))):
        outer_loop_valid = False

    for j in xrange(1, n+1):
        if (diag[j-1] <= 0.0):
            outer_loop_valid = False
     
    
    if (outer_loop_valid):
        iflag[0] = 1
        fcn (n, x, fvec, iflag)
        nfev[0] = 1
        if (iflag[0] < 0):
            return
        fnorm = vector_denorm (n, fvec)
     

    iter = 1 
    ncsuc = 0 
    ncfail = 0 
    nslow1 = 0 
    nslow2 = 0 
    jeval = False 

    while (outer_loop_valid):
        jeval = True

        if (iopt == 1):
            jac(n, x, fvec, fjac, ldfjac, iflag)
            njev[1] = njev[1] + 1
        else:
            iflag[0] = 2
            dfdjc1 (fcn, n, x, fvec, fjac, ldfjac, iflag, ml, mu,
                    epsfcn, wa1, wa2)
            nfev[0] = nfev[0] + min(ml+mu+1, n)
         
        if (iflag[0] < 0):
            outer_loop_valid = False 
            break 
         
        ipvt = [0]
        pivot = False
        
        dqrfrac (n, n, fjac, ldfjac, pivot, ipvt, 1, wa1, wa2, wa3);

        if (iter == 1):
            if (mode == 1):
                for j in xrange(1, n+1):
                    diag[j-1] = wa2[j-1]
                    if (wa2[j-1] == 0.0):
                        diag[j-1] = 1.0
                 
            else:
                for j in xrange(1, n+1):
                    wa3[j-1] = diag[j-1]*x[j-1]
                 
                xnorm = vector_denorm(n, wa3) 
                delta = factor * xnorm 
                if (delta == 0):
                    delta = factor 
                 
             
        
        for i in xrange(1, n+1):
            qtf[i-1] = fvec[i-1] 
         

        for j in xrange(1, n+1):
            if (fjac[j-1][j-1] != 0.0):
                sum_ = 0.0;
                for i in xrange(j, n+1):
                    sum_ += fjac[i-1][j-1]*qtf[i-1]
                 
                temp = -sum_/fjac[j-1][j-1]
                for i in xrange(j, n+1):
                    qtf[i-1] = qtf[i-1] + fjac[i-1][j-1]*temp
                 
             
         
        
        sing = [False]
        for j in xrange(1, n+1):
            l = j
            jm1 = j - 1
            if (jm1 >= 1):
                for i in xrange(1, jm1+1):
                    r[l-1] = fjac[i-1][j-1]
                    l = l + n - i
                 
             
            r[l-1] = wa1[j-1]
            if (wa1[j-1] == 0.0):
                sing[0] = False
             
         
        dqform (n, n, fjac, ldfjac, wa1)

        if (mode == 1):
            for j in xrange(1, n+1):
                diag[j-1] = max(diag[j-1], wa2[j-1])
             
         

        inner_loop = True

        while (inner_loop):

            if (nprint > 0):
                iflag[0] = 0
                if (((iter-1)%nprint) == 0):
                    fcn (n, x, fvec, iflag)
                
                if (iflag[0] < 0):
                    return
             

            ddoglg (n, r, lr, diag, qtf, delta, wa1, wa2, wa3);

            for j in xrange(1, n+1):
                wa1[j-1] = -wa1[j-1]
                
                if bounds is None or bounds[j-1] is None:
                    wa2[j-1] = x[j-1] + wa1[j-1]
                
                elif bounds[j-1] is not None:
                    #We are asked to enforce bounds
                    
                    #Calculate the new_x
                    new_x = x[j-1] + wa1[j-1]
                    
                    lower_bound, upper_bound = bounds[j-1]
                    
                    #Check if the lower bound is violated
                    #if new_x is less than lower_bound ...
                    #change new_x so that it is halfway between the old x
                    #and the lower bound
                    if lower_bound is not None and new_x < lower_bound:
                        #How far did the new_x move?
                        distance = math.sqrt(abs(x[j-1] - new_x))
                        
                        new_x = (lower_bound + x[j-1])/distance
                    
                    #if new_x is bigger than the upper_bound..
                    #change new_x so that it is halfway between the old x
                    #the upper bound
                    elif upper_bound is not None and new_x > upper_bound:
                        distance = math.sqrt(abs(new_x - x[j-1]))
                        
                        new_x = (upper_bound + x[j-1])/distance
                    
                    #Each variable can only violate one variable bound at a 
                    #time. So no problems here
                    
                    #Done with everything
                    wa2[j-1] = new_x
                    
                #End of bounds enforcement
                wa3[j-1] = diag[j-1]*wa1[j-1]
             

            pnorm = vector_denorm(n, wa3)

            if (iter == 1):
                delta = min(delta, pnorm)
             

            iflag[0] = 1
            fcn (n, wa2, wa4, iflag)
            nfev[0] = nfev[0] + 1
            if (iflag[0] < 0):
                return;
            
            fnorm1 = vector_denorm(n, wa4)

            actred = -1.0
            if (fnorm1 < fnorm):
                actred = 1.0 - (fnorm1/fnorm)*(fnorm1/fnorm);
            

            l = 1;
            for i in xrange(1, n+1):
                sum_ = 0.0;
                for j in xrange(i, n+1):
                    #if (l <= len(r)):
                    if (l <= lr):
                        sum_ += r[l-1]*wa1[j-1]
                        l = l + 1
                     
                wa3[i-1] = qtf[i-1] + sum_
             
            temp = vector_denorm (n, wa3)
            prered = 0.0
            if (temp < fnorm):
                prered = 1.0 - (temp/fnorm)*(temp/fnorm);
             
            ratio = 0.0
            if (prered > 0.0):
                ratio = actred/prered;
             

            if (ratio < 0.1):
                ncsuc = 0;
                ncfail = ncfail + 1;
                delta = 0.5*delta;
                
            else:
                ncfail = 0;
                ncsuc = ncsuc + 1;
                if (ratio > 0.5 or ncsuc > 1):
                    delta = max(delta, pnorm/0.5)
                 
                if (abs(ratio - 1.0) <= 0.1):
                    delta = pnorm/0.5;
                 
             

            if (ratio >= 0.0001):

                for j in xrange(1, n+1):
                    x[j-1] = wa2[j-1]
                    wa2[j-1] = diag[j-1] * x[j-1]
                    fvec[j-1] = wa4[j-1]
                 
                xnorm = vector_denorm (n, wa2)
                fnorm = fnorm1
                iter = iter + 1
             

            nslow1 = nslow1 + 1 
            if (actred >= 0.001):
                nslow1 = 0;
             
            if (jeval):
                nslow2 = nslow2 + 1;
             
            if (actred >= 0.1):
                nslow2 = 0;
             

            if (delta <= xtol*xnorm or fnorm == 0.0):
                info[0] = 1;
             
            if (info[0] != 0):
                outer_loop_valid = False;
                break;
             

            if (nfev[0] >= maxfev):
                info[0] = 2;
             
            if (0.1*max(0.1*delta, pnorm) <= EPSMCH*xnorm):
                info[0] = 3;
             
            if (nslow2 == 5):
                info[0] = 4;
             
            if (nslow1 == 10):
                info[0] = 5;
             
            if (info[0] != 0):
                outer_loop_valid = True;
                break;
             
            if (ncfail == 2):
                break;
             
            
            for j in xrange(1, n+1):
                sum_ = 0.0
                for i in xrange(1, n+1):
                    sum_ += fjac[i-1][j-1]*wa4[i-1]
                 
                wa2[j-1] = (sum_ - wa3[j-1])/pnorm
                wa1[j-1] = diag[j-1]*((diag[j-1]*wa1[j-1]))/pnorm
                if (ratio >= 0.0001):
                    qtf[j-1] = sum_
                 
             
            d1updt (n, n, r, lr, wa1, wa2, wa3, sing)
            d1mpyq_m (n, n, fjac, ldfjac, wa2, wa3)
            d1mpyq (1, n, qtf, 1, wa2, wa3)

            jeval = False
         
     
    if (iflag[0] < 0):
        info[0] = iflag[0]
     
    iflag[0] = 0
    
    if (nprint > 0 ):
        fcn (n, x, fvec, iflag)
     
 

def dnsqe (fcn, jac, iopt, n, x, fvec,  tol, nprint, info, bounds=None):
    info[0] = 0;
    
    if (iopt < 1 or iopt > 2 or n <= 0 or tol < 0.0):
        return;
     
    maxfev = 100*(n + 1)
    if (iopt == 2):
        maxfev = 2*maxfev;
     
    xtol = tol
    ml = n - 1
    mu =  n - 1
    epsfcn = 0.0
    mode = 2
    factor = 100.0;
    qtf  = [0.0]*n
    diag = [0.0]*n
    wa1 = [0.0]*n
    wa2 = [0.0]*n
    wa3 = [0.0]*n
    wa4 = [0.0]*n
    fjac = [[0.0]*n for i in range(0, n)]
    nfev = [0]
    njev = [0] 
    for j in xrange(1, n+1):
        diag[j-1] = 1.0
    
    ldfjac = n
    lr = int(math.floor((n*(n + 1))/2))

    r = [0.0]*lr
    
    dnsq (fcn, jac, iopt, n, x, fvec, fjac, ldfjac, xtol,
          maxfev, ml, mu, epsfcn, diag, mode, factor, nprint,
          info, nfev, njev, r, lr, qtf, wa1, wa2, wa3, wa4, bounds);

    if (info[0] == 5):
        info[0] = 4

    return


def dnsqe_nice(fcn, jac, x, tol=math.sqrt(EPSMCH), bounds=None):
    n = len(x)
    f = [0.0]*n
    info = [0]
    iopt = 2
    if jac != None:
        iopt = 1
    nprint = 0
    dnsqe(fcn, jac, iopt, n, x, f,  tol, nprint, info, bounds)
    info_map = {0: "Improper input parameters.\nPlease report this error.",
                1: "Algorithm estimates that the relative error between x and the solution is at most TOL",
                2: "Too many iterations required for solution.\nTry another initial guess.",
                3: "The tolerance is too small.\nNo further improvement in the approximate solution X is possible.\nYour equations are likely badly scaled and very tight\n",
                4: "Iteration is not making good progress.\nThe solver is likely stuck in a variable space valley.\nTry another initial guess.\n"}
    
    fnorm = vector_denorm(n, f)
    return x, f, fnorm, info[0], info_map[info[0]]


def fcn_test1(n, x, fvec, iflag):
    for k in xrange(1, n+1):
        temp = (3.0 - 2.0*x[k-1])*x[k-1]
        temp1 = 0
        if (k != 1):
            temp1 = x[(k-1)-1]
        temp2 = 0
        if (k != n):
            temp2 = x[(k+1)-1]
        fvec[k-1] = temp - temp1 - 2*temp2 + 1.0
     
    return


def test1(n):
    iopt = 2
    x = [-1.0]*n
    fvec = [0.0]*n
    tol = math.sqrt(EPSMCH)
    nprint = 0
    info = [0]
    dnsqe (fcn_test1, None, iopt, n, x, fvec,  tol, nprint, info)    
    return x

def dqrslv (n, r, ldr, ipvt, diag, qtb, x, sigma, wa):

    for j in frange(1,n):
        for i in frange(1, n):
            r[i-1][j-1] = r[j-1][i-1]
         
        x[j-1] = r[j-1][j-1]
        wa[j-1] = qtb[j-1]
     

    for j in frange(1, n):
        l = ipvt[j-1]
        if (diag[l-1] != 0.0):
            for k in frange(j, n):
                sigma[k-1] = 0.0;

            sigma[j-1] = diag[l-1]
            qtbpj = 0.0

            for k in frange(j, n):

                if (sigma[k-1] != 0.0):
                    
                    if (abs(r[k-1][k-1]) < abs(sigma[k-1])):
                        cotan = r[k-1][k-1]/sigma[k-1]
                        sin = 0.5/math.sqrt(0.25 + 0.25*cotan*cotan)
                        cos = sin*cotan
                     
                    else:
                        tan = sigma[k-1]/r[k-1][k-1]
                        cos = 0.5/math.sqrt(0.25 + 0.25*tan*tan)
                        sin = cos*tan
                     

                    r[k-1][k-1] = cos*r[k-1][k-1] + sin*sigma[k-1]
                    temp = cos*wa[k-1] + sin*qtbpj
                    qtbpj = -sin*wa[k-1] + cos*qtbpj
                    wa[k-1] = temp

                    kp1 = k + 1
                    if (n <= kp1):
                        for i in frange(kp1, n):
                            if (i-1) < len(r) and (k-1)< len(r[0]) :
                                temp = cos*r[i-1][k-1] + sin*sigma[i-1]
                                sigma[i-1] = -sin*r[i-1][k-1] + cos*sigma[i-1]
                                r[i-1][k-1] = temp
                         
                     
                 
             
         
        sigma[j-1] = r[j-1][j-1]
        r[j-1][j-1] = x[j-1]
     

    nsing = n
    for j in frange(1, n):
        if (sigma[j-1] == 0.0 and nsing == n):
            nsing = j - 1
         
        if (nsing < n):
            wa[j-1] = 0.0
         
     
    
    if (nsing >= 1):
        for k in frange(1, nsing):
            j = nsing - k + 1
            sum_ = 0.0
            jp1 = j + 1
            if (nsing >= jp1):
                for i in frange(jp1, nsing):
                    sum_ += r[i-1][j-1]*wa[i-1]
                 
             
            wa[j-1] = (wa[j-1] - sum_)/sigma[j-1]
         
    
    for j in frange(1, n):
        l = ipvt[j-1]
        x[l-1] = wa[j-1]

    return


def dmpar (n, r, ldr, ipvt, diag, qtb, delta, par, x, sigma, wa1, wa2):    
    dwarf = DWARF    
    parc = 0.0
    nsing = n
    
    for j in frange(1, n):
        wa1[j-1] = qtb[j-1]
        if (r[j-1][j-1] == 0 and nsing == n):
            nsing = j - 1
        
        if (nsing < n):
            wa1[j-1] = 0.0

    if (nsing >= 1):
        for k in frange(1, nsing):
            j = nsing - k + 1
            wa1[j-1] = wa1[j-1]/r[j-1][j-1]
            temp = wa1[j-1]
            jm1 = j - 1
            if (jm1 >= 1):
                for i in frange(1, jm1):
                    wa1[i-1] = wa1[i-1] - r[i-1][j-1]*temp

    for j in frange(1, n):
        l = ipvt[j-1]
        x[l-1] = wa1[j-1]
    
    iter = 0
    for j in frange(1, n):
        wa2[j-1] = diag[j-1]*x[j-1]
    
    dxnorm = vector_denorm(n, wa2)
    fp = dxnorm - delta
    
    if (fp > 0.1*delta):

        parl = 0.0;
        if (nsing >= n):
            for j in frange(1, n):
                l = ipvt[j-1]
                wa1[j-1] = diag[l-1]*(wa2[l-1]/dxnorm)
            
            for j in frange(1, n):
                sum_ = 0.0
                jm1 = j - 1
                if (jm1 >= 1):
                    for i in frange(1, jm1):
                        sum_ += r[i-1][j-1]*wa1[i-1]
                    
                 
                wa1[j-1] = (wa1[j-1] - sum_)/r[j-1][j-1]
            
            temp = vector_denorm(n, wa1)
            parl = ((fp/delta)/temp)/temp
        
        for j in frange(1, n):
            sum_ = 0.0
            for i in frange(1, j):
                sum_ += r[i-1][j-1]*qtb[i-1]
            
            l = ipvt[j-1]
            wa1[j-1] = sum_/diag[l-1]
        
        
        gnorm = vector_denorm(n, wa1)
        paru = gnorm/delta
        if (paru == 0.0):
            paru = dwarf/min(delta, 0.1)
        

        par[0] = max(par[0], parl)
        par[0] = min(par[0], paru)
        if (paru == 0.0):
            par[0] = gnorm/dxnorm
        

        while (True):

            iter = iter + 1
            if (par[0] == 0.0):
                par[0] = max(DWARF, 0.001*paru)
             
            temp = math.sqrt(par[0])
            for j in frange(1, n):
                wa1[j-1] = temp*diag[j-1]
             
            dqrslv(n, r, ldr, ipvt, wa1, qtb, x, sigma, wa2)

            for j in frange(1, n):
                wa2[j-1] = diag[j-1]*x[j-1]
             
            dxnorm = vector_denorm(n, wa2)
            temp = fp
            fp = dxnorm - delta

            if (abs(fp) <= 0.1*delta or (parl == 0 and fp <= temp and temp <0) or iter == 10):
                break
             

            for j in frange(1, n):
                l = ipvt[j-1]
                wa1[j-1] = diag[l-1]*(wa2[l-1]/dxnorm)
             
            for j in frange(1, n):
                wa1[j-1] = wa1[j-1]/sigma[j-1]
                temp = wa1[j-1]
                jp1 = j + 1
                if (n >= jp1):
                    for i in frange(jp1, n):
                        wa1[i-1] = wa1[i-1] - r[i-1][j-1]*temp;
                    

            if (fp > 0):
                parl = max(parl, par[0])
            
            if (fp < 0):
                paru = min(paru, par[0])
            
            par[0] = max(parl, par[0] + parc);
    
    if (iter == 0):
        par[0] = 0.0

    return



def dwupdt (n, r, ldr, w, b, alpha, cos, sin):

    for j in frange(1, n):
        rowj = w[j-1]
        jm1 = j - 1

        if (jm1 < 1):
            for i in frange(1, jm1):
                temp = cos[i-1]*r[i-1][j-1] + sin[i-1]*rowj
                rowj = -sin[i-1]*r[i-1][j-1] + cos[i-1]*rowj
                r[i-1][j-1] = temp
            

        cos[j-1] = 1.0
        sin[j-1] = 0.0

        if (rowj != 0.0):

            if (abs(r[j-1][j-1]) < abs(rowj)):
                cotan = r[j-1][j-1]/rowj
                sin[j-1] = 0.5/math.sqrt(0.25 + 0.25*(cotan*cotan))
                cos[j-1] = sin[j-1]*cotan
            
            else:
                tan = rowj/r[j-1][j-1]
                cos[j-1] = 0.5/Math.sqrt(0.25 + 0.25*(tan*tan))
                sin[j-1] = cos[j-1]*tan

            r[j-1][j-1] = cos[j-1]*r[j-1][j-1] + sin[j-1]*rowj
            temp = cos[j-1]*b[j-1] + sin[j-1]*alpha
            alpha = -sin[j-1]*b[j-1] + cos[j-1]*alpha
            b[j-1] = temp
        
    

def dfdjc3 (fcn, m, n,  x, fvec, fjac, ldfjac, iflag, epsfcn, wa):
    eps = math.sqrt(max(epsfcn, EPSMCH))
    iflag[0] = 1
    for j in frange(1, n):
        temp = x[j-1]
        h = eps*abs(temp)
        x[j-1] = temp + h
        fcn (iflag, m, n, x, wa, fjac, ldfjac)
        if (iflag[0] < 0):
            return;
        x[j-1] = temp
        for i in frange(1, m):
            fjac[i-1][j-1] = (wa[i-1] - fvec[i-1])/h
        
    


def dnls1 (fcn, jac, iopt, m, n, x, fvec, fjac, ldfjac,
           ftol, xtol, gtol, maxfev, epsfcn, diag, mode,
           factor, nprint, info, nfev, njev, ipvt, qtf,
           wa1, wa2, wa3, wa4):
    
    ijunk = None;
    outer_loop_valid = True
    iflag = [0]
    par = [0.0]
    passthrough1 = False
    info  = [0] 
    iflag = [0]
    nfev = [0]
    njev = [0]
    delta = 0.0
    temp = 0.0
    xnorm = 0.0
    
    if (iopt < 1 or iopt > 3):
        outer_loop_valid = false
    
    if (n <= 0 or m < n):
        outer_loop_valid = false
    
    if (ldfjac < n or ftol < 0.0 or xtol < 0.0 or gtol < 0.0):
        outer_loop_valid = false
    
    if (maxfev <= 0 or factor < 0):
        outer_loop_valid = false
    
    if (iopt < 3 and ldfjac < m):
        outer_loop_valid = false
    
    if (mode == 2):
        for j in frange(1, n):
            if (diag[j-1] <= 0.0):
                outer_loop_valid = false
        
    iflag[0] = 1
    ijunk = 1
    fcn(iflag, m, n, x, fvec, ijunk, ijunk)
    nfev[0] = 1

    if (iflag[0] < 0):
        outer_loop_valid = false
    

    fnorm = vector_denorm(m, fvec)
    par[0] = 0.0
    iter = 1

    while (outer_loop_valid):

        if (nprint > 0):
            iflag[0] = 0
            if (((nprint)%(iter-1)) == 0):
                fcn(iflag, m, n, x, fvec, fjac, ijunk)
            
            if (iflag[0] < 1):
                outer_loop_valid = false
                break
            
        
        
        if (iopt == 2):
            iflag[0] = 2
            fcn(iflag, m, n, x, fvec, ijunk, ijunk)
            njev[0] = njev[0] + 1
            #Skip derivative check
            passthrough1 = True
        

        # Code approximates the Jacobian
        if (iopt == 1 or passthrough1):
            if (not passthrough1):
                iflag[0] = 1
                dfdjc3 (fcn, m, n, x, fvec, fjac, ldfjac, iflag, epsfcn, wa4)
                nfev[0] = nfev[0] + n
           
            if (iflag[0] < 0):
                outer_loop_valid = False
                break
           
            dqrfrac (m, n, fjac, ldfjac, True, ipvt, n, wa1, wa2, wa3)
            for i in frange(1, m):
                wa4[i-1] = fvec[i-1]
           
            for j in frange(1, n):
                if (fjac[j-1][j-1] != 0.0):
                    sum = 0.0
                    for i in frange(j, m):
                        sum += fjac[i-1][j-1]*wa4[i-1]
                    temp = -sum/fjac[j-1][j-1]
                    for i in frange(j, m):
                        wa4[i-1] = wa4[i-1] + fjac[i-1][j-1]*temp
                    
                fjac[j-1][j-1] = wa1[j-1]
                qtf[j-1] = wa4[j-1]
           
        

        #Calculate the Jacobian Matrix
        if (iopt == 3):
            for j in frange(1, n):
                qtf[j-1] = 0.0
                for i in frange(1, n):
                    fjac[i-1][j-1] = 0.0
                
            
            for i in frange(1, m):
                nrow = i
                iflag[0] = 3
                fcn(iflag, m, n, x, fvec, wa3, nrow)
                if (iflag[0] < 0):
                    outer_loop_valid = False
                    break;
                
                #Skip the derivative check
                temp = fvec[i-1]
                dwupdt (n, fjac, ldfjac, wa3, qtf, temp, wa1, wa2)
            
            njev[0] = njev[0] + 1

            sing = False
            for j in frange(1, n):
                if (fjac[j-1][j-1] == 0):
                    sing = True
                
                ipvt[j-1] = j
                wa2[j-1] = column_denorm(j, fjac, 1, j)
            
            if(sing):
                dqrfrac (n, n, fjac, ldfjac, true, ipvt, n, wa1, wa2, wa3);
                for j in frange(1, n):
                    if (fjac[j-1][j-1] != 0.0):
                        sum = 0.0;
                        for i in frange(j, n):
                            sum += fjac[i-1][j-1]*qtf[i-1]
                        
                        temp = -sum/fjac[j-1][j-1]
                        for i in frange(j, n):
                            qtf[i-1] = qtf[i-1] + fjac[i-1][j-1]*temp
                        
                    
                    fjac[j-1][j-1] = wa1[j-1]
                
            
        

        #Line 560
        if (iter == 1):
            if (mode != 2):
                for j in frange(1, n):
                    diag[j-1] = wa2[j-1]
                    if (wa2[j-1] == 0.0):
                        diag[j-1] = 1.0;
                    
                
            
            for j in frange(1, n):
                wa3[j-1] = diag[j-1]*x[j-1]
            
            xnorm = vector_denorm(n, wa3)
            delta = factor*xnorm
            if (delta == 0.0):
                delta = factor
            
        
        gnorm = 0.0
        if (fnorm != 0.0):
            for j in frange(1, n):
                l = ipvt[j-1]
                if (wa2[l-1] != 0.0):
                    sum = 0.0
                    for i in frange(1, j):
                        sum += fjac[i-1][j-1]*(qtf[i-1]/fnorm)
                    
                    gnorm = max(gnorm, abs(sum/wa2[l-1]))
                
            
        
        if (gnorm <= gtol):
            info[0] = 4
        
        if (info[0] != 0):
            outer_loop_valid = False
            break
        

        if (mode != 2):
            for j in frange(1, n):
                diag[j-1] = max(diag[j-1], wa2[j-1])
            
        
        while(True):

            dmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, par, wa1, wa2, wa3, wa4);

            for j in frange(1, n):
                wa1[j-1] = -wa1[j-1]
                wa2[j-1] = x[j-1] + wa1[j-1]
                wa3[j-1] = diag[j-1]*wa1[j-1]
            
            pnorm = vector_denorm(n, wa3)

            if (iter == 1):
                delta = min(delta, pnorm)

            iflag[0] = 1
            fcn(iflag, m, n, wa2, wa4, fjac, ijunk)
            nfev[0] = nfev[0] + 1
            if (iflag[0] < 0):
                outer_loop_valid = false
                break
            
            fnorm1 = vector_denorm(m, wa4)

            actred = -1.0
            if (0.1*fnorm1 < fnorm):
                actred = 1.0 - (fnorm1/fnorm)*(fnorm1/fnorm)

            for j in frange(1, n):
                wa3[j-1] = 0.0
                l = ipvt[j-1]
                temp = wa1[l-1]
                for i in frange(1, j):
                    wa3[i-1] = wa3[i-1] + fjac[i-1][j-1]*temp
                
            
            temp1 = vector_denorm(n, wa3)/fnorm
            temp2 = math.sqrt(par[0])*pnorm/fnorm
            prered = temp1*temp1 + temp2*temp2/0.5
            dirder = -(temp1*temp1 + temp2*temp2)
            
            ratio = 0.0;
            if (prered != 0.0):
                ratio = actred/prered

            if (ratio <= 0.25):
                if (actred >= 0.0):
                    temp = 0.5
                if (actred < 0.0):
                    temp = 0.5*dirder/(dirder + 0.5*actred)
                if (0.1*fnorm1 >= fnorm or temp < 0.1):
                    temp = 0.1
                delta = temp*min(delta, pnorm/0.1)
                par[0] = par[0]/temp
                
            else:
                if (par[0] == 0.0 or ratio >= 0.75):
                    delta = pnorm/0.5;
                    par[0] = 0.5*par[0];
            
            if (ratio >= 0.0001):
                
                for j in frange(1, n):
                    x[j-1] = wa2[j-1]
                    wa2[j-1] = diag[j-1]*x[j-1]
                
                for i in frange(1, m):
                    fvec[i-1] = wa4[i-1]

                xnorm = vector_denorm(n, wa2);
                fnorm = fnorm1;
                iter = iter + 1;
            
            if (abs(actred) <= ftol and prered <= ftol and 0.5*ratio <= 1.0):
                info[0] = 1
            
            if (delta <= xtol*xnorm):
                info[0] = 2
            
            if (abs(actred) <= ftol and prered <= ftol and 0.5*ratio <= 1 and info[0] == 2):
                info[0] = 3
            
            if (info[0] != 0):
                outer_loop_valid = False
                break

            if (nfev[0] >= maxfev):
                info[0] = 5
            
            if (abs(actred) <= EPSMCH and prered <= EPSMCH and 0.5*ratio <= 1.0):
                info[0] = 6
            
            if (delta <= EPSMCH*xnorm):
                info[0] = 7

            if (gnorm <= EPSMCH):
                info[0] = 8
            
            if (info[0] != 0):
                break
            
            if (ratio < 0.0001):
                #do nothing complete the inner-loop
                pass
            else:
                break
            
        
        if (not outer_loop_valid):
            break
        

def dnls1e(fcn, jac, iopt, m, n, x, fvec, tol, nprint, info):
    ftol = tol;
    xtol = tol;
    gtol = 0.0;
    epsfcn = 0.0;
    mode = 1;
    maxfev = 100*(n + 1);
    ldfjac = 0
    
    diag = [0.0]*n
    nfev = [0]
    njev = [0]
    ipvt = [0]*n
    qtf = [0.0]*n
    wa1 = [0.0]*n
    wa2 = [0.0]*n
    wa3 = [0.0]*n
    wa4 = [0.0]*m
    factor = 100.0
    
    fjac = [[0.0]*n for i in range(0, m)]
    
    if (iopt == 1):
        maxfev = maxfev*2
        fjac = [[0.0]*n for i in range(0, m)]
        ldfjac = m
        
    else:
        #Is this supposed to be nxn?
        fjac = [[0.0]*n for i in range(0, m)]
        ldfjac = n;

    dnls1(fcn, jac, iopt, m, n, x, fvec, fjac, ldfjac,
          ftol, xtol, gtol, maxfev, epsfcn, diag, mode,
          factor, nprint, info, nfev, njev, ipvt, qtf,
          wa1, wa2, wa3, wa4)

def f_nle(iflag, m, n, x, f, fjac, ldfjac):
    id_011 = 0.7 #B
    id_012 = 1.7 #A
    id_003 = 0.2 #z1
    id_006 = 0.8 #z2
    id_008 = 88.32 #t
    id_000 = x[0] #x1
    id_001 = x[1] #alpha
    id_002 = x[2] #k1
    id_004 = x[3] #x2
    id_005 = x[4] #k2
    id_007 = x[5] #p1
    id_009 = x[6] #p2
    id_010 = x[7] #gamma2
    id_013 = x[8] #gamma1
    id_014 = x[9] #y1
    id_015 = x[10] #y2
    f[0] = (id_000*((1+(id_001*((id_002-(1))))))-(id_003)) - (0) + (1-math.exp(id_001))
    f[1] = (id_004*((1+(id_001*((id_005-(1))))))-(id_006)) - (0)
    f[2] = (id_000+(id_004)) - (1)
    f[3] = (id_007) - (math.pow(10,(7.62231-(1417.9/((191.15+(id_008)))))))
    f[4] = (id_009) - (math.pow(10,(8.10765-(1750.29/((235+(id_008)))))))
    f[5] = (id_010) - (math.pow(10,(id_011*(id_000)*(id_000)/((math.pow((id_000+(id_011*(id_004)/(id_012))),2))))))
    f[6] = (id_013) - (math.pow(10,(id_012*(id_004)*(id_004)/((math.pow((id_012*(id_000)/(id_011)+(id_004)),2))))))
    f[7] = (id_002) - (id_013*(id_007)/(760))
    f[8] = (id_005) - (id_010*(id_009)/(760))
    f[9] = (id_014) - (id_002*(id_000))
    f[10] = (id_015) - (id_005*(id_004))
    #f[11] = (1-math.exp(id_001))
    return True


def fcn(iflag, m, n, x, fvec, fjac, ldfjac):
    y = [1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1,
         3.5e-1, 3.9e-1, 3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1,
         1.34, 2.1, 4.39]

    for i in frange(1, m):
        tmp1 = i
        tmp2 = 16.0 - i
        tmp3 = tmp1
        if (i > 8):
            tmp3 = tmp2
        fvec[i-1] = y[i-1] - (x[1-1] + tmp1/(x[2-1]*tmp2 + x[3-1]*tmp3))

def test2():    
    n = 11
    m = 11
    x = [1.0]*n
    x[0] = 0.1
    x[1] = 1.0
    x[2] = 0.9
    fvec = [0.0]*m
    tol = math.sqrt(EPSMCH)
    iopt = 1
    nprint = 0
    info = [0]
    dnls1e(f_nle, None, iopt, m, n, x, fvec, tol, nprint, info)
    print 'Results:'
    print '\n'.join([str(v) for v in x])
    



if __name__ == '__main__':
    #import cProfile
    #cProfile.run("test1(50)")
    #print test1(9)
    test2()
