def gras(tabin=None,coltot=None,rowtot=None,iter_in=None):


    import numpy as np
# # Generated with SMOP  0.41
# from libsmop import *
# # gras.m

    
# @function

    # #ras
# balance tabin to column total (coltot) and row total (rowtot)
# input tabin, coltot, rowtot
    
    # tabin= initial estimate of table to be balanced
# coltot= column totals to be reached of the table (must be same dimensions of table
# rowtot= row totals to be reached of the table (must be same dimensions of table
# iter_in = number of interations (optional)
    
    #optional check for dimensionality:
# if length(coltot)~=size(tabin,2)
#     disp('warning sizes dont match')
# elseif length(rowtot)~=size(tabin,1)
#     disp('warning sizes dont match')
# end
    
    # check that the totals match, otherwise it will never converge!
    if abs(abs(sum(rowtot)) - abs(sum(coltot))) > 1e-07:
        print('row and col totals do not match')
        abs(sum(rowtot))
        abs(sum(coltot))
    
# % split the input table into a positive and negative 
    postab = (tabin>=0)*tabin
    negtab = -1*(tabin<0)*tabin
    
    # just set up some dimension variables
    tabdim1=postab.shape[0]
    tabdim2=postab.shape[1]

    rdim=len(rowtot)
    
    sdim=len(coltot)
    
    # if the maximum number of iterations is externally defined (nargin=4), use
# it, otherwise use 100 iterations
    # if nargin == 4:
    #     MAXITER=iter_in
    # else:
    MAXITER=100

    
    print('MAXITER: ' + str(MAXITER))
    #initialise the row and column scaling vectors (r1 and s1) to unity
    # if nargin > 2:
        # r1=np.ones((rdim,MAXITER),'float')
   
    r1=np.ones((rdim,MAXITER+1),'float')

    # s1=ones(MAXITER,length(coltot))
    s1=np.ones((MAXITER+1,sdim),'float')

    for k in range(1,MAXITER):
        tmp=scalcer(postab,negtab,r1[0:rdim,k],coltot,sdim,tabdim1)
        s1[k,0:sdim] = tmp.transpose()
        r1[0:rdim,k+1] = rcalcer(postab,negtab,s1[k,0:sdim],rowtot,rdim,tabdim2)
 
        if (sum(abs(r1[:,k+1] - r1[:,k])) < 1e-08):
            if k > 1:
                if (sum(abs(s1[k,:]- s1[k-1,:])) < 1e-08):
                    print('threshold reached: ' + k)
                    break
    
        if k == MAXITER:
            print('max runs reached: ' + str(k))
            break
    
    
    r2=(abs(r1) < 1e-10) + r1
    s2=(abs(s1) < 1e-10) + s1
    # r2=min(r2,10000000000.0)
    # s2=min(s2,10000000000.0)

    # % postab = diag(r1(:,k+1))*[postab(1:size(r1,1),1:size(s1,2))*diag(s1(k,:)),postab(1:size(r1,1),size(s1,2)+1:size(postab,2))];
    postab1 = np.matmul(postab,np.diag(s1[k,:]))
    postab2 = np.matmul(np.diag(r1[:,k+1]),postab1)

    # % negtab = inv(diag(r2(:,k+1)))*[negtab(1:size(r1,1),1:size(s1,2))*inv(diag(s2(k,:))),negtab(1:size(r1,1),size(s1,2)+1:size(negtab,2))];
    # negtab1 = np.matmul(negtab,np.linalg.inv(np.diag(max(1e-5,s2[k,:]))))
    negtab1 = np.multiply(negtab,1/s2[k,:])
    
    # negtab2 = np.matmul(np.linalg.inv(np.diag(max(1e-5,r2[:,k+1]))),negtab1)
    negtab2  = np.multiply(np.vstack(1/r2[:,k+1]),negtab1)

    tabout = postab2 - negtab2
    
    np.var*22
    
    return tabout


def rcalcer(p=None,n=None,s=None,u=None,dim1=None,dim2=None,*args,**kwargs):
#Calculate the scaling factor on rows
#p is positive matrix
#n is negative matrix
#s is scalar
#dim1 is no of rows of p/n
#dim2 is no of cols of p/n
   import numpy as np
   
   r=np.ones((dim1),'float')

   for i in range(1,dim1):
        pscale=rpcalcer(p,s,i,dim2)
        nscale=rncalcer(n,s,i,dim2)
        if pscale == 0:
            if nscale == 0:
                r[i]=1
            else:                
                r[i]=(u[i] + np.sqrt((u[i]) ** 2 + (4*nscale))) / (2)
        else:
            r[i]=(u[i] + np.sqrt((u[i]) ** 2 + 4*pscale*nscale)) / (2*pscale)
        return r
    
def rpcalcer(p=None,s=None,i=None,dim2=None,*args,**kwargs):
    pofr=0
    for j in range(1,dim2):
        if j >len(s):
            pofr=pofr + p[i,j]
        else:
            pofr=pofr + (p[i,j]*s[j])
    return pofr
    
def rncalcer(n=None,s=None,i=None,dim2=None):

    nofr=0
    for j in range(1,dim2):
        if j > len(s):
            nofr=nofr + n[i,j]
        else:
            nofr=nofr + n[i,j] / (s[j] + 1e-10)
    return nofr

    
def scalcer(p=None,n=None,r=None,v=None,dim1=None,dim2=None,*args,**kwargs):

    import numpy as np
    
#Calculate the scaling factor on columns
#p is positive matrix
#n is negative matrix
#s is scalar
#dim1 is no of rows of p/n
#dim2 is no of cols of p/n
    
    s=np.ones((dim1),'float')
    for j in range(1,dim1):
        pscale=spcalcer(p,r,j,dim2)
        nscale=sncalcer(n,r,j,dim2)
        if pscale == 0:
            if nscale == 0:
                s[j]=1
            else:
                s[j]=(v[j] + np.sqrt((v[j]) ** 2 + 4*nscale)) / (2)
        else:
            s[j]=(v[j] + np.sqrt((v[j]) ** 2 + 4*pscale*nscale)) / (2*pscale)
    
    return s
    
       
def spcalcer(p=None,r=None,j=None,dim2=None):
    
    pofr=0
    for i in range(1,dim2):
        if i > len(r):
            pofr=pofr + p[i,j]
        else:
            pofr=pofr + (p[i,j]*r[i])
    return pofr
    
def sncalcer(n=None,r=None,j=None,dim2=None):
    nofr=0
    for i in range(1,dim2):
        if i > len(r):
            nofr=nofr + n[i,j]
        else:
            nofr=nofr + n[i,j] / (r[i] + 1e-10)
    return nofr

    