function [clusters] =creat_clusters_fixk_kmeans(data,H)
[n,~]=size(data);
clusters=zeros(n,H);
maxk=floor(sqrt(n));
k=min(maxk,50);
for i=1:H
    clusters(:,i)=kmeans(data,k,'emptyaction','singleton');
end

end