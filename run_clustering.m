function [acc,NMI]=run_clustering(Z,gnd)
%,ACC,NMI
K=max(gnd);
%% 
Z = ( abs(Z) + abs(Z') ) / 2 ;
Z = processC(Z,0.9);  %  0.5:0.1:0.9

% refining Z
[U,S,V] = svd(Z);
S = diag(S);
r = min(4*K+1,sum(S>1e-3*S(1)));
S = S(1:r);
U = U(:,1:r)*diag(sqrt(S));
U = normr(U);
Z = U*U';Z=abs(Z);
L = Z.^4;

idx = clu_ncut(L,K);  % turnable
NMI=100*MutualInfo(gnd,idx);
acc = 100*compacc(idx,gnd);





