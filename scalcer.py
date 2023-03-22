# Generated with SMOP  0.41
from libsmop import *
# scalcer.m

    
@function
def scalcer(p=None,n=None,r=None,v=None,dim1=None,dim2=None,*args,**kwargs):
    varargin = scalcer.varargin
    nargin = scalcer.nargin

    # #ras
#Calculate the scaling factor on columns
#p is positive matrix
#n is negative matrix
#s is scalar
#dim1 is no of rows of p/n
#dim2 is no of cols of p/n
    
    s=ones(dim1,1)
# scalcer.m:10
    for j in arange(1,dim1).reshape(-1):
        pscale=pcalcer(p,r,j,dim2)
# scalcer.m:13
        nscale=ncalcer(n,r,j,dim2)
# scalcer.m:14
        if pscale == 0:
            if nscale == 0:
                s[j]=1
# scalcer.m:17
            else:
                s[j]=(v(j) + sqrt((v(j)) ** 2 + dot(4,nscale))) / (2)
# scalcer.m:19
        else:
            s[j]=(v(j) + sqrt((v(j)) ** 2 + dot(dot(4,pscale),nscale))) / (dot(2,pscale))
# scalcer.m:22
            #       s(j) = (exp(1)*v(j) + sqrt((exp(1)*v(j))^2+4*pscale*nscale))/(2*pscale);
    
    
    
    
    
    
@function
def pcalcer(p=None,r=None,j=None,dim2=None,*args,**kwargs):
    varargin = pcalcer.varargin
    nargin = pcalcer.nargin

    pofr=0
# scalcer.m:32
    for i in arange(1,dim2).reshape(-1):
        if i > size(r,1):
            pofr=pofr + p(i,j)
# scalcer.m:35
        else:
            pofr=pofr + dot(p(i,j),r(i))
# scalcer.m:37
    
    
@function
def ncalcer(n=None,r=None,j=None,dim2=None,*args,**kwargs):
    varargin = ncalcer.varargin
    nargin = ncalcer.nargin

    nofr=0
# scalcer.m:44
    for i in arange(1,dim2).reshape(-1):
        #     if r(i) ~=0
        if i > size(r,1):
            nofr=nofr + n(i,j)
# scalcer.m:48
        else:
            nofr=nofr + n(i,j) / (r(i) + 1e-10)
# scalcer.m:50
        #     end
    