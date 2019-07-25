function [RendUnd,DriftCorr] = DriftCorrection(RendUnd,RiskPremium,rfree,d,sigmaUnd,sigmaCurr,rho,WrkDaysR,WrkDays,Drift)

% Function that modifies drift inside a simulated distribution

% Internal Variable
    NumSim = size(RendUnd,1);
    NumAsset = size(RendUnd,2);
    NumSteps = size(RendUnd,3);
    
% Memory Space
    DriftCorr = zeros(1,NumAsset,NumSteps);
    
% Drift Correction
    meanRendUnd = mean(RendUnd,1);
    RendUnd = RendUnd-repmat(meanRendUnd,NumSim,1,1);
    dt = diff(WrkDays);
    dt = [WrkDays(1);dt];
    if strcmp(Drift,'RF')   % Risk Free correction
        for i = 1:NumSteps
            for j = 1:NumAsset
                %DriftCorr(1,j,i) = ((rfree(i,j)-d(i,j))./WrkDaysR-0.5*sigmaUnd(i,j).^2-rho(j).*sigmaUnd(j).*sigmaCurr(j)).*dt(i);
                DriftCorr(1,j,i) = ((rfree(i,j)-d(i,j))./WrkDaysR-0.5*sigmaUnd(i,j).^2-rho(i,j).*sigmaUnd(i,j).*sigmaCurr(i,j)).*dt(i);
            end
        end
        DriftCorr = cumsum(DriftCorr,3);
    elseif strcmp(Drift,'HR')   % Historical Return correction
        for i = 1:NumSteps
            for j = 1:NumAsset
                %DriftCorr(1,j,i) = (RiskPremium(i,j)-0.5*sigmaUnd(i,j).^2).*dt(i);
                DriftCorr(1,j,i) = RiskPremium(i,j).*dt(i);
            end
        end
        DriftCorr = cumsum(DriftCorr,3);
    elseif strcmp(Drift,'ER')   % Expected Return correction
        for i = 1:NumSteps
            for j = 1:NumAsset
                DriftCorr(1,j,i) = ((rfree(i)+RiskPremium(i)-d(j))./WrkDaysR-0.5*sigmaUnd(j).^2-rho(j).*sigmaUnd(j).*sigmaCurr(j)).*dt(i);
            end
        end
        DriftCorr = cumsum(DriftCorr,3);
    end
    RendUnd = RendUnd+repmat(DriftCorr,NumSim,1,1);