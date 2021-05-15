%% label reconstruction
% this function is used to change the label of SVM datasets
function y=label_reconstruction(y, name)
    if strcmp(name, 'mnist.scale')
        y(y~=0)=1;
        y(y==0)=-1;        
    elseif strcmp(name, 'cifar10')
        y(y~=0)=1;
        y(y==0)=-1;        
    elseif strcmp(name, 'SVHN.scale')
        y(y==1)=-1;        
        y(y>1)=1;
    elseif strcmp(name, 'covtype.libsvm.binary.scale')
        y(y==1)=-1;
        y(y==2)=1;       
    elseif strcmp(name, 'news20.scale')
        y(y==1)=-1;
        y(y>1)=1;       
    elseif strcmp(name, 'sector.scale')
        y(y~=1)=1;
        y(y==1)=-1;
    end

end