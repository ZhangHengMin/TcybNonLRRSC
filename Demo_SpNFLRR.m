clear all;
close all; 
clc;

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
paras.La=[5 6 7 8 9 10]; % [5 6 7 8 9 10]turnable 
paras.ranks= 100;%[50, 80, 100, 130, 150];   % turnable 
method={'SpNF_LRR1','SpNF_LRR23','SpNF_LRR12'};
%% select method

for ll = 1 : length(method)
    fun = method{ll};
    disp([ 'method = ' num2str(fun)]);
    %% our method
    for i = 1:length(paras.Mu)
        disp([ 'mu = ' num2str(paras.Mu(i))]); 
        for k = 1:length(paras.La)
             
            %
            for j = 1:length(paras.ranks)
               
                tic;
                switch fun 
                    case 'SpNF_LRR1'
                        [Z,iter,~] = SpNF_LRR1(X,paras.La(k),paras.ranks(j),paras.Mu(i));

                    case 'SpNF_LRR23'
                        [Z,iter,~] = SpNF_LRR23(X,paras.La(k),paras.ranks(j),paras.Mu(i));

                    case 'SpNF_LRR12'
                        [Z,iter,~] = SpNF_LRR12(X,paras.La(k),paras.ranks(j),paras.Mu(i));
                end
                Time= toc;
                % Results
                [acc,NMI]=run_clustering(Z,gnd); % Improved results 

                disp([ ' lambda_value = ' num2str(paras.La(k)),' acc_result= ' num2str(acc,'%.2f') ,  ' NMI_result= ' num2str(NMI,'%.2f'), ...
                    ' time_cost = ' num2str(Time,'%.2f'),' iter_cost = ' num2str(iter) ]);
            end

        end
    end
end
 
