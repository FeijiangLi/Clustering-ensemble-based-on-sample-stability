function [result] =DA_quadratic_simple_link(clusters,k,da)
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

dataed=predata(da);
Y1=dataed;


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

figure();
for i=1:detial_k
    aa=liu((find(result==i)),1);
%     plot(Y1(aa,1),Y1(aa,2),'.','markersize',30);
    plot3(Y1(aa,1),Y1(aa,2),Y1(aa,3),'.','markersize',30);
%     axis square;
    hold all
end
set(gca,'XLim',[0 1]);
set(gca,'YLim',[0 1]);
set(gca,'FontSize',26);
% saveas(gcf,'4.fig')
% saveas(gcf,'4.eps','psc2')



[newliu,newliudata]=de_result3(locat,liu,de,sim,detial_k,data,Y1);
data(newliu,end+1)=newliudata(:,end);
result=data(:,end);
detial_k=length(unique(result));



figure();
for i=1:detial_k
    aa= result==i;
%     plot(Y1(aa,1),Y1(aa,2),'.','markersize',30);
    plot3(Y1(aa,1),Y1(aa,2),Y1(aa,3),'.','markersize',30);
%     axis square;
    hold all     
end
set(gca,'XLim',[0 1]);
set(gca,'YLim',[0 1]);
set(gca,'FontSize',26);
% saveas(gcf,'5.fig')
% saveas(gcf,'5.eps','psc2')


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

figure();
for i=1:k
    aa= result==i;
%     plot(Y1(aa,1),Y1(aa,2),'.','markersize',30);
    plot3(Y1(aa,1),Y1(aa,2),Y1(aa,3),'.','markersize',30);
%     axis square
    hold all     
end
set(gca,'XLim',[0 1]);
set(gca,'YLim',[0 1]);
set(gca,'FontSize',26);
% saveas(gcf,'6.fig')
% saveas(gcf,'6.eps','psc2')


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
end

function [binvat] =step1(vat)
th=graythresh(vat);
binvat=vat>th;
end