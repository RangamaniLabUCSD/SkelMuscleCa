%% start parameter optimization
VOnly = false;
withMito = true;
if VOnly
    startFromPrev = true;
    % highSensIdxOld = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
    highSensIdx = [1,3,4,5,8,9,11,13,14,16,18,19,22,23,24,25,26,28,33,40,62,77,79,80,81,82,99];
    highSensIdxOld = [1,3,5,8,9,11,13,14,16,18,22,23,24,25,26,28,33,40,46,76,77,79,80,81,82];
    % highSensIdx = [1,4,5,6,9,11,13,14,16,18,19,22,23,24,25,26,30,33,40,47,76,77,79,80,81,82,98];
    if startFromPrev
        % load Data/pVec_VOnly.mat pVec
        VOnlyStruct = load('pSol_fullWithMito/PSO_05-Jul-2026_withMito_Vonly.mat');%load('Data/pVec_VOnly.mat');
        pVec = VOnlyStruct.pSol;
        pTest = ones(1,106);
        pTest(highSensIdxOld) = pVec;
        initVec = pTest(highSensIdx);
    else
        initVec = ones(1,length(highSensIdx));
    end
    lb = 0.5*ones(length(highSensIdx),1); 
    ub = 2*ones(length(highSensIdx),1);
else
    highSensIdx = 1:106;
    VOnlyIdx = [1,3,4,5,8,9,11,13,14,16,18,19,22,23,24,25,26,28,33,40,62,77,79,80,81,82,99];
    VOnlyIdxOld = [1,3,5,8,9,11,13,14,16,18,22,23,24,25,26,28,33,40,46,76,77,79,80,81,82];
    % VOnlyIdx = [1,4,5,6,9,11,13,14,16,18,19,22,23,24,25,26,30,33,40,47,76,77,79,80,81,82,98];
    % VOnlyIdxOld = [1,3,4,5,6,8,9,11,13,14,16,18,19,22,23,24,25,26,28,30,33,40,76,77,79,80,81,82];
    VOnlyIdxCur = VOnlyIdxOld;
    % load Data/pVec_VOnlyNew.mat pVec
    VOnlyStruct = load('pSol_fullWithMito/PSO_05-Jul-2026_withMito_Vonly.mat');%load('Data/pVec_VOnly.mat');
    pVec_VOnly = VOnlyStruct.pSol;
    notInclIdx = setdiff(VOnlyIdxOld,VOnlyIdx);
    for idx = notInclIdx
        pVec_VOnly(VOnlyIdxOld==idx) = 1; % don't allow these non sens to start from diff baseline
    end
    % VOnlyStruct = load('pSol_fullWithMito/pBest_withMito_Vonly_Jul14.mat');%load('Data/pVec_VOnly.mat');
    % pVec_VOnly = VOnlyStruct.pVec;
    varyFactor = 2;
    startFromPrev = true;
    lb = (1/varyFactor)*ones(length(highSensIdx),1);
    ub = varyFactor*ones(length(highSensIdx),1);
    if startFromPrev
        % load pSol_allParam.mat pSol
        load pSol_fullWithMito/PSO_12-Jul-2026_withMito.mat pSol
        initVec = pSol;
    else
        initVec = ones(1,length(highSensIdx));
    end
    for i = 1:length(lb)
        if any(VOnlyIdxCur == highSensIdx(i))
            lb(i) = pVec_VOnly(VOnlyIdxCur == highSensIdx(i))*(1/varyFactor);
            ub(i) = pVec_VOnly(VOnlyIdxCur == highSensIdx(i))*varyFactor;
        end
    end
end
% Setting bounds for parameters

load p0Struct.mat p0Struct
p0 = p0Struct.data;
if withMito
    p0 = [p0; 1; 1; 1];
end
for i = 1:length(highSensIdx)
    if any(highSensIdx(i)==76:83) % voltage params
        if VOnly
            lower_lim = p0(highSensIdx(i))-10;
            upper_lim = p0(highSensIdx(i))+10;
        else
            if any(VOnlyIdxCur == highSensIdx(i))
                prefactor = pVec_VOnly(VOnlyIdxCur == highSensIdx(i));
            else
                prefactor = 1;
            end
            % prefactor = 1;
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
% if any(initVec > ub) || any(initVec < lb)
%     fprintf('do something')
% end
% Particle Swarm Optimization
[pSol,fval,exitflag] = SkelMuscleCa_paramEst(lb, ub, withMito, initVec);