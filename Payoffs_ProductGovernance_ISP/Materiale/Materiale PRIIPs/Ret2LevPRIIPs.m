    function [Px] = Ret2LevPRIIPs(Meanshift,Rend,Px0,ParallelShift)
        
    % Internal Variable
        NumSim = size(Rend,1);
        NumSteps = size(Rend,3);
        MeanshiftCorr(:,1,:)=repmat(Meanshift',[NumSim,1]);
        Px = repmat(Px0(end,:),[NumSim,1,NumSteps]).*exp(Rend+MeanshiftCorr)-repmat(ParallelShift,[NumSim,1,NumSteps]);
