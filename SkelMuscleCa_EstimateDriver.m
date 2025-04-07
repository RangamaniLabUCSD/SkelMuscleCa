%% start parameter optimization
VOnly = false;
if VOnly
    highSensIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,28,30,32,33,40,76,77,78,79,80,81,82,83,91];
    lb = 0.5*ones(length(highSensIdx),1); 
    ub = 2*ones(length(highSensIdx),1);
else
    highSensIdx = [3,15,18,23,32,35,40,42,43,69,72,77,79,80,81,83,86,90,91];
    VOnlyIdx = [1,3,4,5,6,8,9,10,11,13,14,15,16,17,18,19,22,23,24,25,26,28,30,32,33,40,76,77,78,79,80,81,82,83,91];
    load pVec_VOnly.mat pVec
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
[pSol_LM,fval,exitflag] = SkelMuscleCa_paramEst(lb, ub);
username = getenv('USER');
fname = sprintf('/tscc/lustre/ddn/scratch/%s/psofig_%s.png',username,datetime('today'));
saveas(gcf, fname);
