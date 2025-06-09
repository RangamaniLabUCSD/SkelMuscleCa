function QOI = SkelMuscleCa_SAOutput(param)
% Output Function for Senstivity Analysis
% Output:
%   QOI are the quantities of interst
% Input:
%   param is the list of variable parameters

load Data/p0Struct.mat p0Struct
parameters = p0Struct.data;

totSize = size(param,1);
QOI = zeros(size(param,1),14*3+1);
parpool(50);
parfor (i = 1 : length(param(:,1)))
    pVec = param(i,:)'.*parameters;
    [objVal1, qoiList] = SkelMuscleObj(pVec); % objVal is a scalar, qoiList is a vector of length 14*3
    QOI(i,:) = [qoiList, objVal1]; % 14*3 + 1 in length
    fprintf('Session %d of %d \n',i,totSize);
end
end