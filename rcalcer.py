# Generated with SMOP  0.41
from libsmop import *
# rcalcer.m

    
@function
def rcalcer(p=None,n=None,s=None,u=None,dim1=None,dim2=None,*args,**kwargs):
    varargin = rcalcer.varargin
    nargin = rcalcer.nargin

    # #ras
#Calculate the scaling factor on rows
#p is positive matrix
#n is negative matrix
#s is scalar
#dim1 is no of rows of p/n
#dim2 is no of cols of p/n
    
    r=ones(dim1,1)
# rcalcer.m:10
    for i in arange(1,dim1).reshape(-1):
        pscale=pcalcer(p,s,i,dim2)
# rcalcer.m:13
        nscale=ncalcer(n,s,i,dim2)
# rcalcer.m:14
        if pscale == 0:
            if nscale == 0:
                r[i]=1
# rcalcer.m:17
            else:
                r[i]=(u(i) + sqrt((u(i)) ** 2 + dot(4,nscale))) / (2)
# rcalcer.m:19
        else:
            r[i]=(u(i) + sqrt((u(i)) ** 2 + dot(dot(4,pscale),nscale))) / (dot(2,pscale))
# rcalcer.m:22
            #       r(i) = (exp(1)*u(i) + sqrt((exp(1)*u(i))^2+4*pscale*nscale))/(2*pscale);
    
    
    
    
    
    
@function
def pcalcer(p=None,s=None,i=None,dim2=None,*args,**kwargs):
    varargin = pcalcer.varargin
    nargin = pcalcer.nargin

    pofr=0
# rcalcer.m:32
    for j in arange(1,dim2).reshape(-1):
        if j > size(s,2):
            pofr=pofr + p(i,j)
# rcalcer.m:35
        else:
            pofr=pofr + dot(p(i,j),s(j))
# rcalcer.m:37
    
    
@function
def ncalcer(n=None,s=None,i=None,dim2=None,*args,**kwargs):
    varargin = ncalcer.varargin
    nargin = ncalcer.nargin

    nofr=0
# rcalcer.m:44
    for j in arange(1,dim2).reshape(-1):
        #     if s(j) ~=0
        if j > size(s,2):
            nofr=nofr + n(i,j)
# rcalcer.m:48
        else:
            nofr=nofr + n(i,j) / (s(j) + 1e-10)
# rcalcer.m:50
        #     end
    