% clear all;
% clc

intervel = 100;
save_position = zeros(1000, 2);

for i = 0 : 9
    index = find(gnd == i);
    save_position((i * intervel + 1) : (100 + i * intervel), 1) = index(1:100);
    save_position((i * intervel + 1) : (100 + i * intervel), 2) = i;
end