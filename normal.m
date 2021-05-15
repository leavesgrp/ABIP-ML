clc;
clear;

load('./datasets/ADLITTLE.mat');
% data=struct();
% data.A = sparse(Problem.A); 
% data.b = full(Problem.b);
% data.c = full(Problem.aux.c); 

data = struct('A', sparse(Problem.A), 'b', full(Problem.b), 'c', full(Problem.c));



params = struct();

[x,y,s,info] = abip(data,params);