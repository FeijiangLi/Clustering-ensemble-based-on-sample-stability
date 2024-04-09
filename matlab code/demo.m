load Iris.mat
H=50;
gt=data(:,end);
k=length(unique(gt));
data_feature=data(:,1:end-1);
data_feature=predata(data_feature);
[clusterings] =creat_clusters_fixk_kmeans(data_feature,H);
[cl1] =DA_neighbor_simple_link2(clusterings,k);
[ac1,ARI1,NMI1]=evaluate2(cl1,gt,k)
[cl2] =DA_quadratic_simple_link2(clusterings,k);
[ac2,ARI2,NMI2]=evaluate2(cl2,gt,k)



load chainlink.mat
gt=data(:,end);
k=length(unique(gt));
data_feature=data(:,1:end-1);
data_feature=predata(data_feature);
[clusterings] =creat_clusters_fixk_kmeans(data_feature,H);
[cl3] =DA_neighbor_simple_link(clusterings,k,data_feature);
[ac3,ARI3,NMI3]=evaluate2(cl3,gt,k)
[cl4] =DA_quadratic_simple_link(clusterings,k,data_feature);
[ac4,ARI4,NMI4]=evaluate2(cl4,gt,k)
