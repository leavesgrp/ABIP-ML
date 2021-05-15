%% L1 SVM experiment
clear;
clc; 
close all;
diary off

addpath('abip');

Probname = {'australian'};
nprob = length(Probname);
Problist = [1:nprob];

alpha=1.7;
adapt=[5];
normalize=1;
scalar=[1];


for i=1:length(alpha)
    for j=1:length(adapt)
        params = struct('Problem', 'L1_SVM' );

        for di = 1:length(Problist) 
            probID = Problist(di);
            name = Probname{probID};
            [y,X]=libsvmread(['./datasets/',name]);
            y=label_reconstruction(y, name);

            for k=1:length(scalar)   
                data.X=X;
                data.y=y;
                data.scalar=scalar(k);

%                 folder=['./L1_SVM_results/alpha',mat2str(alpha),'/',name,'/'];
%                 if exist(folder)==0
%                     mkdir(folder); 
%                 end

%                 diary(['./L1_SVM_results/alpha',mat2str(alpha),'/',name,'/alpha',mat2str(alpha(i)),'adapt',mat2str(adapt(j)),'scalar',mat2str(data.scalar),...
%                                   'normalize',mat2str(normalize),'.txt']);
%                 diary on
                [x,~,~,~] = abip(data,params);
%                 diary off
%                 clc;
            end
%             clear data X y
%             clc;
        end
    end
end