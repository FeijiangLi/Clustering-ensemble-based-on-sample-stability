function [newliu,newliudata]=de_result3(locat,liu,de,sim,k,data,Y1)
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

%     relocati=sum(liu,2);
%     relocati=relocati./sum(l);
   
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
   [newliu,newliudata]=de_result3(newliu,newliudata,newde,sim,k,data,Y1);
end


end

function [binvat] =step1(vat)
th=graythresh(vat);
binvat=vat>th;
end