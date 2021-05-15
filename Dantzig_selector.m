%% Dantzig Selector experiment
clear;
clc; 
close all;
diary off

addpath('abip');

Probname = {'pyrim'};
nprob = length(Probname);

Problist = [1:nprob];



alpha=1.8;
adapt=20;
method=1;
normalize=0;
scalar = [0.1];


for i=1:length(alpha)
    for j=1:length(adapt)
        
        params = struct('Problem', 'Dantzig_selector' );

        for di = 1:length(Problist) 
            probID = Problist(di);
            name = Probname{probID};
            [y,X]=libsvmread(['./datasets/',name]);
            for k=1:length(scalar)
                data.X=X;
                data.y=y;
                data.scalar=scalar(k);

%                 folder=['./Dantzig_selector_results/',name,'/'];
%                 if exist(folder)==0
%                     mkdir(folder); 
%                 end

%                 diary([folder,'alpha',mat2str(alpha(i)),'adapt',mat2str(adapt(j)),'scalar',mat2str(data.scalar),...
%                                   'normalize',mat2str(normalize),'METHOD',mat2str(method),'.txt']);
%                 diary on
                [x,~,~,info] = abip(data,params);
%                 diary off
%                 clc
            end    
%             clear data X y
%             clc;
        end
    end
end
