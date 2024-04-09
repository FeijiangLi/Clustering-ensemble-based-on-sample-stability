function [result] =DA_quadratic_simple_link2(clusters,k)
n=size(clusters,1);
data=(1:1:n)';
[sim]=simnumber(clusters);
th=graythresh(sim);
sim1=sim>=th;
sim1=sim1.*(1-th);

sim2=sim<th;
sim2=sim2.*(th);

sim12=sim1+sim2;

SH=((sim-th)./sim12).^2;
he=sum(SH)./n;
he=(he-min(he))./(max(he)-min(he));
th=graythresh(he);
locat=find(he>th);

A=1:1:n;
de=setdiff(A,locat);
liu=data(locat,:);
newsim=sim(locat,locat);
newsim_dist=1-newsim;
newsim_dist=newsim_dist-diag(diag(newsim_dist));
newsim_dist = squareform(newsim_dist);
Z=linkage(newsim_dist,'single');
I=inconsistent(Z,4);
[~,find_k]=max(I(:,end));
detial_k=length(locat)+1-find_k(end);
if detial_k<k
    detial_k=k;  
end
result=cluster(Z,'maxclust',detial_k);
liu=[liu result];

[newliu,newliudata]=de_result3(locat,liu,de,sim,detial_k,data);
data(newliu,end+1)=newliudata(:,end);
result=data(:,end);
detial_k=length(unique(result));

if detial_k~=k
    simij=zeros(1,(detial_k-1)*(detial_k-2)/2+(detial_k-1));
    for i=1:detial_k
        locat_i=result==i;
        for j=i+1:detial_k
            locat_j=result==j;
            sim_ij=sim(locat_i,locat_j);
            ij=sum(sum(sim_ij))/(sum(locat_i)*sum(locat_j)); 
            simij((detial_k-1)*(i-1)-((i-1)*(i-2)/2)+j-i)=ij;
        end
    end
    simij=simij./max(simij);
    dis=1-simij;
    Z = linkage(dis,'single');
    clu=cluster(Z,k);
    CC=zeros(n,1);
    for i=1:detial_k
        biaoji=clu(i);
        l= result==i;
        CC(l,:)=biaoji;
    end
    result=CC;
end


function [sim] = simnumber(clusters) 
[N,M] = size(clusters); 
newE = zeros(N,M);
ucl = unique(clusters(:,1)); 
if (max(clusters(:,1) ~= length(ucl)))
    for j = 1:length(ucl)
        newE(clusters(:,1) == ucl(j),1) = j; 
    end
end

for i = 2:M
    ucl = unique(clusters(:,i)); 
    prevCl = length(unique(newE(:,[1:i-1])));
    for j = 1:length(ucl)
        newE(clusters(:,i) == ucl(j),i) = prevCl + j; 
    end
end
no_allcl = max(max(newE));

[n,m]=size(newE);
max_E=max(newE);
min_E=min(newE);
relabel=zeros(n,no_allcl);
for i=1:m
    cl=newE(:,i);
    for j=min_E(i):max_E(i)
        relocat=cl==j;
        relabel(relocat,j)=1;
    end
end
sim=relabel*relabel'./M;



function [newliu,newliudata]=de_result3(locat,liu,de,sim,k,data)
newliudata=liu;
deliu=sim(de,locat);
n=length(de);
relocat=zeros(n,k);
for i=1:k
    l= newliudata(:,end)==i;
    liu=deliu(:,l);
    max_liu=max(max(liu));
    min_liu=min(min(liu));
    if max_liu==min_liu
        liul=liu;
        liul(:)=1;
    else
        liul=(liu-min(min(liu)))./(max(max(liu))-min(min(liu)));
    end
    binliul=step1(liul);
    binliul_value=liu.*binliul;
    relocati=sum(binliul_value,2);
    ni=sum(binliul,2)+1;
    relocati=relocati./ni;
    relocat(:,i)=relocati;
end

if sum(sum(relocat))==0
    sumzeros=ones(size(relocat,1),1);
    result=sumzeros.*(k+1);
    k=k+1;
else
    maxrelocat=max(relocat,[],2);
    maxrelocat=maxrelocat./max(maxrelocat);
    binsumrelocat=step1(maxrelocat);
    relocat=relocat(binsumrelocat,:);
    sumzeros=binsumrelocat==1;
    [~,result]=max(relocat,[],2);
end

newliu=de(sumzeros);
newde=setdiff(de,newliu);
newdeliudata=data(newliu,:);
newliu=[locat newliu];
newdeliudata(:,end+1)=result;
newliudata=[newliudata;newdeliudata];

if ~isempty(newde)
   [newliu,newliudata]=de_result3(newliu,newliudata,newde,sim,k,data);
end

function [binvat] =step1(vat)
th=graythresh(vat);
binvat=vat>th;
