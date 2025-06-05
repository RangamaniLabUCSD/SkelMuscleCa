%% start parameter optimization
VOnly = false;
if VOnly
    highSensIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
    lb = 0.5*ones(length(highSensIdx),1); 
    ub = 2*ones(length(highSensIdx),1);
else
    highSensIdx = 1:106;
    VOnlyIdx = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
    load pVec_VOnlyNoBib.mat pVec
    lb = 0.5*ones(length(highSensIdx),1);
    ub = 2*ones(length(highSensIdx),1);
    for i = 1:length(lb)
        if any(VOnlyIdx == highSensIdx(i))
            lb(i) = pVec(VOnlyIdx == highSensIdx(i))/2;
            ub(i) = pVec(VOnlyIdx == highSensIdx(i))*2;
        end
    end
end
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