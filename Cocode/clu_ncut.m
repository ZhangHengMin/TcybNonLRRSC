function [idx] = clu_ncut(L,K)
% this routine groups the data X into K subspaces by NCut
% inputs:
%       L -- an N*N affinity matrix, N is the number of data points
%       K -- the number of subpaces (i.e., clusters)
L = (abs(L) + abs(L'))/2;
%imshow(L);
D = diag(1./sqrt(sum(L,2)));
L = D*L*D;
% 
[U,S,V] = svd(L);
V = U(:,1:K);
V = D*V; % laplace graph
rand('twister',5489); %% fixed 
idx = kmeans(V,K,'emptyaction','singleton','replicates',100,'display','off');
idx = idx';