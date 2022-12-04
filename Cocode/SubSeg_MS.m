function acc = SubSeg_MS(SegMethod, X, gnd, lambda)

switch SegMethod
    
    case 'lrr'  % special code
        [L,E] = lrr(X,lambda);
        
    case 'lrra'  % special code
        [L,E] = solve_lrr(X,lambda);
        
        %--------------------Schatten-p minimization norm based LRR--------------------------
    case 'spdualNulrr'
        [L,E] = lrr(X,lambda);
        
    case 'spdual12lrr'
        [L,E] = nlrr_dualS12(X,lambda);
        
    case 'spdual23lrr'
        [L,E] = nlrr_dualS23(X,lambda);
        
        %--------------------Schatten-p decomposition norm based LRR--------------------------
 
    case 'spdualNudlrr'
        [L,E] = lrr(X,lambda);
        
    case 'spdual12dlrr'
        [L,E] = nlrr_dualS12(X,lambda);
        
    case 'spdual23dlrr'
        [L,E] = nlrr_dualS23(X,lambda);
end

% % 下面的代码对MS有降性能作用 对FC有升性能作用
% for i = 1 : size(L,2)
%    L(:,i) = L(:,i) / max(abs(L(:,i))) ;
% end   % normalized coefficient matrix

nCluster = length(unique(gnd));
idx = clu_ncut(L,nCluster);  % turnable

% seg_results
acc = compacc(idx,gnd);  % different evaluation metric
