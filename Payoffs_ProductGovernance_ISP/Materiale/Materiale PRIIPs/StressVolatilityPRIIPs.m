function [vStress,Stat] = StressVolatilityPRIIPs(Data,LogRend,T,Freq,UnderlyingType)    

% Parameter    
    if strcmp(Freq,'D')
        if T <= 1
            VolaWindow = 21;
        elseif T > 1
            VolaWindow = 63;
        end
    elseif strcmp(Freq,'W')
        if T <= 1
            VolaWindow = 8;
        elseif T > 1
            VolaWindow = 16;
        end
    elseif strcmp(Freq,'M')
        if T <= 1
            VolaWindow = 6;
        elseif T > 1
            VolaWindow = 12;
        end
    end
    if T <= 1
        VolaPerc = 99;
    elseif T > 1
        VolaPerc = 90;
    end

% Internal Variable
    NumAsset = size(LogRend,2);
    NumDataRet = size(LogRend,1);
    

    if strcmp(UnderlyingType{1,1},'Rate')
        
    vStress = zeros(1,NumAsset);
    
    for i=1:NumAsset

    Data1 = Data(Data(:,i*2-1)>0,i*2);
    LogRendFull = log(Data1(2:end,1)./Data1(1:end-1,1));
    
% Memory Space    
    NumDataRet = size(LogRendFull,1);
    StepRoll = NumDataRet-VolaWindow;
    Stat = zeros(StepRoll+1,NumAsset);   

% Volatility Distribution
    for t = 0:StepRoll
    % Window Selection
        DataLog = LogRendFull(1+t:VolaWindow+t,:);
    % Volatility   
        v = std(DataLog,1);
        Stat(t+1,i) = v;
    end
    Stat = sort(Stat);
    Percentile = (size(Stat,1)+1)*VolaPerc/100;
    vStress(i) = Stat(floor(Percentile),i)+(Percentile-floor(Percentile))*(Stat(ceil(Percentile),i)-Stat(floor(Percentile),i));
    end
    
    % Stat solo per grafico finale
    
    % Memory Space    
    NumDataRet = size(LogRend,1);
    StepRoll = NumDataRet-VolaWindow;
    Stat = zeros(StepRoll+1,NumAsset);    
        
% Volatility Distribution
    for t = 0:StepRoll
    % Window Selection
        DataLog = LogRend(1+t:VolaWindow+t,:);
    % Volatility   
        v = std(DataLog,1);
        Stat(t+1,:) = v;
    end

    else

    % Memory Space    
    NumDataRet = size(LogRend,1);
    StepRoll = NumDataRet-VolaWindow;
    Stat = zeros(StepRoll+1,NumAsset);    
        
% Volatility Distribution
    for t = 0:StepRoll
    % Window Selection
        DataLog = LogRend(1+t:VolaWindow+t,:);
    % Volatility   
        v = std(DataLog,1);
        Stat(t+1,:) = v;
    end
    Stat = sort(Stat);
     Percentile = (size(Stat,1)+1)*VolaPerc/100;
    if Percentile>size(Stat,1)
        Percentile=size(Stat,1);
    end
    vStress = Stat(floor(Percentile),:)+(Percentile-floor(Percentile))*(Stat(ceil(Percentile),:)-Stat(floor(Percentile),:));
    end
    
%vStress = prctile(Stat,VolaPerc);
