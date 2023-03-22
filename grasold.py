def gras(tabin, coltot,rowtot,iter_in):

# function [tabout,r1,s1] = gras(tabin, coltot,rowtot,iter_in) %#eml
# % #ras
# % balance tabin to column total (coltot) and row total (rowtot)
# % input tabin, coltot, rowtot

# % tabin= initial estimate of table to be balanced
# % coltot= column totals to be reached of the table (must be same dimensions of table
# % rowtot= row totals to be reached of the table (must be same dimensions of table
# % iter_in = number of interations (optional)


# %optional check for dimensionality:
# % if length(coltot)~=size(tabin,2)
# %     disp('warning sizes dont match')
# % elseif length(rowtot)~=size(tabin,1)
# %     disp('warning sizes dont match')
# % end


# % check that the totals match, otherwise it will never converge!
    if (abs(abs(sum(rowtot))-abs(sum(coltot)))> 1e-7) :
        print('row and col totals do not match')
        abs(sum(rowtot))
        abs(sum(coltot))


# % split the input table into a positive and negative 
    postab = (tabin>=0)*tabin
    negtab = -(tabin<0)*tabin


# % just set up some dimension variables
    tabdim1=size(postab,1)
    tabdim2=size(postab,2)
    # if nargin>2,rdim=length(rowtot);end
    sdim=length(coltot)

# % if the maximum number of iterations is externally defined (nargin=4), use
# % it, otherwise use 100 iterations
    if nargin==4:
        MAXITER=iter_in
    else:
        MAXITER=100

    print(['MAXITER: ',num2str(MAXITER)])

# initialise the row and column scaling vectors (r1 and s1) to unity
    if nargin>2:
        r1 = ones(length(rowtot),MAXITER)
        
    s1 = ones(MAXITER,length(coltot));


    for k in range(1,MAXITER):
        s1(k,:) = scalcer(postab,negtab,r1(:,k),coltot,sdim,tabdim1)
        r1(:,k+1) = rcalcer(postab,negtab,s1(k,:),rowtot,rdim,tabdim2)
 
        if (sum(abs(r1(:,k+1)-r1(:,k)))<1e-8)   
            if k>1
                if (sum(abs(s1(k,:)-s1(k-1,:)))<1e-8) 
                    fprintf('threshold reached: ')
                    fprintf(num2str(k))
                    fprintf('\n')
                    break
                
    if k==MAXITER
        fprintf('max runs reached: ')
        fprintf(num2str(k))
        fprintf('\n')
# %     r1(:,k+1)
# %     r1(:,k)
    
    
    r2 = (abs(r1)<1e-10)+r1
    s2 = (abs(s1)<1e-10)+s1
    r2= min(r2,1e10)
    s2= min(s2,1e10)

    # % postab = diag(r1(:,k+1))*[postab(1:size(r1,1),1:size(s1,2))*diag(s1(k,:)),postab(1:size(r1,1),size(s1,2)+1:size(postab,2))];
    postab1 = [postab(:,1:size(s1,2))*diag(s1(k,:)),postab(1:size(postab,1),size(s1,2)+1:size(postab,2))];
    postab2 = [diag(r1(:,k+1))*postab1(1:size(r1,1),:);postab1(size(r1,1)+1:size(postab1,1),1:size(postab1,2))];

    # % negtab = inv(diag(r2(:,k+1)))*[negtab(1:size(r1,1),1:size(s1,2))*inv(diag(s2(k,:))),negtab(1:size(r1,1),size(s1,2)+1:size(negtab,2))];
    negtab1 = [negtab(:,1:size(s1,2))*inv(diag(max(1e-5,s2(k,:)))),negtab(1:size(negtab,1),size(s1,2)+1:size(negtab,2))];
    negtab2 = [inv(diag(max(1e-5,r2(:,k+1))))*negtab1(1:size(r1,1),:);negtab1(size(r1,1)+1:size(negtab1,1),1:size(negtab,2))];

    tabout = postab2 - negtab2;
    
    
def rcalcer(p,n,s,u,dim1,dim2) %#eml
    % #ras
    %Calculate the scaling factor on rows
    %p is positive matrix
    %n is negative matrix
    %s is scalar
    %dim1 is no of rows of p/n
    %dim2 is no of cols of p/n

    r=ones(dim1,1);

    for i=1:dim1
        pscale = pcalcer(p,s,i,dim2);
        nscale = ncalcer(n,s,i,dim2);
        if pscale==0
            if nscale==0
                r(i) = 1;
            else
                r(i) = (u(i) + sqrt((u(i))^2+4*nscale))/(2);
            end
        else
            r(i) = (u(i) + sqrt((u(i))^2+4*pscale*nscale))/(2*pscale);
    %       r(i) = (exp(1)*u(i) + sqrt((exp(1)*u(i))^2+4*pscale*nscale))/(2*pscale);
        end
    end
        
        
        
        
    function pofr = pcalcer(p,s,i,dim2)

    pofr = 0;
    for j =1:dim2
        if j>size(s,2)
            pofr = pofr +p(i,j);
        else
            pofr = pofr + p(i,j)*s(j);
        end 
    end


    function nofr = ncalcer(n,s,i,dim2)

    nofr = 0;
    for j =1:dim2
        %     if s(j) ~=0
        if j>size(s,2)
            nofr = nofr+n(i,j);
        else
            nofr = nofr + n(i,j)/(s(j)+1e-10);
        end 
        %     end
    end