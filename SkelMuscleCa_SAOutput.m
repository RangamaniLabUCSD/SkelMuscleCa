function QOI = SkelMuscleCa_SAOutput(param)
% Output Function for Senstivity Analysis
% Output:
%   QOI are the quantities of interst
% Input:
%   param is the list of variable parameters

input = importdata('InputParam1.xlsx'); % Load default parameters
parameters = input.data; % Parameter values
% load PSO_25-Apr-2024.mat yinit


totSize = size(param,1);
QOI = zeros(size(param,1),14*4+2);
parpool(50);
parfor (i = 1 : length(param(:,1)))
    pVec = param(i,:)'.*parameters;
    [objVal1, qoiList] = SkelMuscleObj(pVec); % objVal is a scalar, qoiList is a vector of length 14*4
    objVal2 = SkelMuscleObj2(pVec);
    QOI(i,:) = [qoiList, objVal1, objVal2]; % 14*4 + 2 in length
    fprintf('Session %d of %d \n',i,totSize);
end
end