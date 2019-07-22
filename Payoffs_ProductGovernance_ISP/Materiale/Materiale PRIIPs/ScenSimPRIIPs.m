function [Meanshift,RendUnd0,RendUnd] = ScenSimPRIIPs(LogRend,Steps,vStress,StressScen,SimType,NumSim,Seed)
% Historical Bootstrapping

% Seed Setting
    if Seed == 'Y'
        reset(RandStream.getGlobalStream)
    end

% Internal Variable
    NumData = size(LogRend,1);
    NumAsset = size(LogRend,2);
    NumSteps = size(Steps,1);
% Memory Space    
    LogRendScen = zeros(NumSim,NumAsset);
    RendUnd = zeros(NumSim,NumAsset,NumSteps);
    Meanshift=zeros(NumData,size(LogRend,2));
    ti = 1; 
% Simulation Model (RTS Annex II. PART I, 5/6)    
    if strcmp(StressScen,'Y')
        LogRend = LogRend-repmat(mean(LogRend,1),NumData,1);
        LogRend = (LogRend./repmat(std(LogRend),NumData,1)).*repmat(vStress,NumData,1);
    end

% Bootstrapping (RTS Annex II. PART I, 22(b))
    if strcmp(SimType{1,1},'PCA')
    % 23.a (v) covariance matrix 
         Meanshift=repmat(mean(LogRend,1),NumData,1);
         LogRend=LogRend- Meanshift;
        varcov = cov(LogRend);
    % 23.a (vi) PCA
        coef = pcacov(varcov);
    % 23.a (vii-viii) NewMat (Select principal component)
          NewMat = coef(:,1:3);
    % 23.a (ix) return projection
        LogRend = LogRend*NewMat;
    % 23.a (x) 
        LogRend = LogRend*NewMat'; 

    end
% 23.b (i) Bootstrapping
    for t = 1:Steps(end)
        pointer = randi([1 size(LogRend,1)],NumSim,1);
        LogRendScen_t = LogRend(pointer,:);
        LogRendScen = LogRendScen+LogRendScen_t;
        if t == Steps(ti)
            RendUnd(:,:,ti) = LogRendScen;
            ti = ti+1;
        end 
    end 
    
   Meanshift=Meanshift(1:NumSteps,:).*repmat(Steps,1,NumAsset); 
% Remove Historical mean
    RendUnd0 = RendUnd-repmat(mean(RendUnd,1),NumSim,1,1);


    