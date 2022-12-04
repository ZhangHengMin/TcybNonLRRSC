clear all;
clc;

%data processing
load('.\0_Data\digital\MNIST_2k2k.mat');

clear testIdx trainIdx;
intervel = 100;
gnd0 = zeros(1000, 2);

for i = 0 : 9
    index = find(gnd == i);
    gnd0((i * intervel + 1) : (100 + i * intervel), 1) = index(1:100);
    gnd0((i * intervel + 1) : (100 + i * intervel), 2) = i+1;
end
    feaA = fea(gnd0(:,1),:);
    gndB = gnd0(:,2);
    
    save ('minst_part', 'feaA', 'gndB');










