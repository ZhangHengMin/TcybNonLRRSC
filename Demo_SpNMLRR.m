clear all;
close all;
 

addpath('.\Dataset');
addpath('.\PROPACK');
addpath('.\Cocode');
addpath('.\SpDLRR');
addpath('.\SpMLRR');
%
dataset='yaleb10';
data=load([dataset,'.mat']);
data=data.obj;
X = data.X;        %m*n
% normalize
for ii = 1 : size(X,2) 
    X(:,ii) = X(:,ii) /norm(X(:,ii)) ;
end
gnd = data.cids;   %1*n

%% para
paras.Mu= [1e-6 1e-4]; % 10^{-8,...,+6}
paras.La=[0.5 1.0 1.5 2.0 3.0]; % turnable
method={'SpNM_LRR1','SpNM_LRR23','SpNM_LRR12'};
%% select method

for ll = 1 : length(method)
    fun = method{ll};
    disp([ 'method = ' num2str(fun)]);
    %% our method
    for i = 1:length(paras.Mu)
        
        disp([ 'mu = ' num2str(paras.Mu(i))]);
        for k = 1:length(paras.La)

            tic;
            switch fun
                case 'SpNM_LRR1'
                    [Z,iter,~] = SpNM_LRR1(X,paras.La(k), paras.Mu(i));

                case 'SpNM_LRR23'
                    [Z,iter,~] = SpNM_LRR23(X,paras.La(k), paras.Mu(i));

                case 'SpNM_LRR12'
                    [Z,iter,~] = SpNM_LRR12(X,paras.La(k), paras.Mu(i));
            end
            Time= toc;
            % Results
            [acc,NMI]=run_clustering(Z,gnd); % Improved results

            disp([ ' lambda_value = ' num2str(paras.La(k)),' acc_result= ' num2str(acc,'%.2f') ,  ' NMI_result= ' num2str(NMI,'%.2f'), ...
                ' time_cost = ' num2str(Time,'%.2f'),' iter_cost = ' num2str(iter) ]);
        end

    end
end

 

