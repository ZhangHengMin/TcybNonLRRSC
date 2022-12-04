function [Z, Iter, E] = SpNF_LRR12(D, lambda, ranks, mu, rho, tol, maxIter)
%
% min_{U1,U2,M1,M2,Z,E}: 1/2*(\|M1\|_{*}+\|M2\|_{*})+lambda*\|E\|_l,
%                                       s.t. AZ + E = D, L = U1U2', U1 = M1£¬U2 = M2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%    D:        Input data matrix of size m*n
%    ranks:    Given rank
%    lambda:   Regularization paramter
%    tol:      Tolerance for stopping criterion
%    maxIter:  Maximum iteration number, usually 500.
%    rho:      Initial value of the parameter
%    mu:     Penalty parameter, it should be small enough
%Output:
%    Z : n*n   Low-rank component
%    E : m*n   Error component
%    errList : Difference

% Default parameters
[ncr(1,1), ncr(1,2)] = size(D);
if nargin < 7
    maxIter = 1000;
end
if nargin < 6
    tol = 1e-4;
end
if nargin < 5
    rho = 1.1;   % turnable 1.05
end
if nargin < 4
    mu = 1e-6;
end
% if nargin < 3
%     ranks = round(1.2*rank_estimation(D));
% end
% if nargin < 2
%     lambda = sqrt(max(ncr));
% end

%%% Initialization
alpha    = [0.5,0.5];
max_beta = 1e10;

for i  = 1:2
    if i == 1
        M{i} = rand(size(D, i), ranks);            % fix
        [U{i}, aa1] = qr(D'*M{i}, 0);
        Y{i} = zeros(size(D, 3-i), ranks);
    else
        M{i} = rand(size(D, 3-i), ranks);
        [U{i}, aa2] = qr(D'*M{i}, 0);
        Y{i} = zeros(size(D, i), ranks);
    end
    %U{i} = eye(size(T, i), ranks); 
    M{i} = Y{i};
end
Z  = U{1}*U{2}';
a1 = norm(Z, 2);   
a2 = norm(Z, Inf)/lambda;
Y2 = Z/max(a1,a2);
Y1 = zeros(ncr(1,1), ncr(1,2));
clear aa a1 a2;

A = D;
AAEye = inv(A'*A + eye(ncr(1,2)));
E  = sparse(ncr(1,1), ncr(1,2));
stop_temp = zeros(2, 1);

%% Main loop
for Iter = 1: maxIter
   %
    beta = 1/mu;
    % Update U1 and U2
    EE   = Z - Y2*beta;
    U{1} = (EE*U{2} + M{1} + Y{1}*beta)*inv(U{2}'*U{2} + eye(ranks));
    U{2} = (EE'*U{1} + M{2} + Y{2}*beta)*inv(eye(ranks) + U{1}'*U{1});
    
    % Update M1 and M2
    for i = 1:2
        [UU, sigma, VV] = svd(U{i} - Y{i}*beta, 'econ');
        sigma = diag(sigma);
        svp = length(find(sigma > alpha(i)*beta));
        if svp>=1
            sigma = sigma(1:svp) - alpha(i)*beta;
        else
            svp = 1;
            sigma = 0;
        end
        M{i} = UU(:,1:svp)*diag(sigma(1:svp))*VV(:,1:svp)';
    end
    
    % Update Z
    X = U{1}*U{2}';
    Z = AAEye*(X + Y2*beta - A'*(E - D + Y1*beta));
    
    % Update E
    AZ = A*Z;
    E_temp = D - AZ - Y1*beta;
    Twothird = 3;
    if Twothird == 1
        E = solve_lf(E_temp, lambda*beta);        % F-norm based operator
    elseif Twothird == 2
        E = solve_l1(E_temp, lambda*beta);        % L1-norm based operator
    else
        E = solve_l1l2(E_temp, lambda*beta);      % L21-norm based operator
    end
    
    % Update Lagrange multipliers
    for i = 1:2
        temp = M{i} - U{i};
        Y{i} = Y{i} + mu*temp;
        stop_temp(i) = norm(temp(:), 'fro');
    end
    Y2 = Y2 + mu*(X - Z);
    leq1 = AZ + E - D;
    Y1 = Y1 + mu*leq1;
    
    %
    stopC = max(max(abs(leq1)));
    if stopC < tol
        break;
    else
        mu = min(mu * rho, max_beta);
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end
function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end

% This Function to Estimate the Rank of the Input Matrix
function d = rank_estimation(X)

[n, m] = size(X);
epsilon = nnz(X)/sqrt(m*n);
mm = min(100, min(m, n));
S0 = lansvd(X, mm, 'L');

S1 = S0(1:end-1)-S0(2:end);
S1_ = S1./mean(S1(end-10:end));
r1 = 0;
lam = 0.05;
while(r1 <= 0)
    for idx = 1:length(S1_)
        cost(idx) = lam*max(S1_(idx:end)) + idx;
    end
    [v2, i2] = min(cost);
    r1 = max(i2-1);
    lam = lam + 0.05;
end
clear cost;

for idx = 1:length(S0)-1
    cost(idx) = (S0(idx+1)+sqrt(idx*epsilon)*S0(1)/epsilon )/S0(idx);
end
[v2, i2] = min(cost);
r2 = max(i2);

d = max([r1 r2]);