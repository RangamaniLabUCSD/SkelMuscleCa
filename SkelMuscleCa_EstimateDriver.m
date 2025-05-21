%% start parameter optimization
VOnly = false;
if VOnly
    % highSensIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,28,30,32,33,40,76,77,78,79,80,81,82,83,91];
    highSensIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
    lb = 0.5*ones(length(highSensIdx),1); 
    ub = 2*ones(length(highSensIdx),1);
else
    highSensIdx = 1:106;%[3,15,18,23,32,35,40,42,43,69,72,77,79,80,81,83,86,90,91];
    % highSensIdx = [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,...
    %     26, 27, 28, 29, 32, 33, 34, 36, 37, 38, 39, 40, 42, 43, 44, 76, 77, 78, 79, 80, 81, 82, 83, 85, 88, 90, 92];
    VOnlyIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
    % VOnlyIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,28,30,32,33,40,76,77,78,79,80,81,82,83,91];
    load pVec_VOnlyNoBib2.mat pVec
    % load PSOpenaltyNoBib_VOnlyRedo_13-May-2025.mat pVec
    % oldIdx = 1:105;
    % load PSOpenaltyEst3_07-Apr-2025.mat pSol
    lb = 0.5*ones(length(highSensIdx),1);
    ub = 2*ones(length(highSensIdx),1);
    for i = 1:length(lb)
        if any(VOnlyIdx == highSensIdx(i))
            lb(i) = pVec(VOnlyIdx == highSensIdx(i))/2;
            ub(i) = pVec(VOnlyIdx == highSensIdx(i))*2;
        end
    end
end
% [1:75,84:105];%[12,31,34,41,44,55:57,59:68,70,71,73,74,75,84,89,93,94,95];
% Setting bounds for parameters

load p0Struct.mat p0Struct
p0 = p0Struct.data;
for i = 1:length(highSensIdx)
    if any(highSensIdx(i)==76:83) % voltage params
        if VOnly
            lower_lim = p0(highSensIdx(i))-10;
            upper_lim = p0(highSensIdx(i))+10;
        else
            if any(VOnlyIdx == highSensIdx(i))
                prefactor = pVec(VOnlyIdx == highSensIdx(i));
            else
                prefactor = 1;
            end
            lower_lim = prefactor*p0(highSensIdx(i))-10;
            upper_lim = prefactor*p0(highSensIdx(i))+10;
        end
        if p0(highSensIdx(i)) < 0
            ub(i) = lower_lim/p0(highSensIdx(i));
            lb(i) = upper_lim/p0(highSensIdx(i));
        else
            lb(i) = lower_lim/p0(highSensIdx(i));
            ub(i) = upper_lim/p0(highSensIdx(i));
        end
    end
end

% Particle Swarm Optimization
[pSol,fval,exitflag] = SkelMuscleCa_paramEst(lb, ub);