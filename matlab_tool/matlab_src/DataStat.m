function [HistMean,HistSigma,HistRho,NumAsset,Pxt0] = DataStat(LogRend,Data,UnderlyingType)

% Internal Variable
    NumAsset = size(Data,2)/2;
    
    if strcmp(UnderlyingType{1,1},'Rate')
    for i=1:NumAsset
    Data1=Data(Data(:,i*2-1)>0,i*2);
    LogRendFull = log(Data1(2:end,1)./Data1(1:end-1,1)); 
    HistSigma(i) = std(LogRendFull,1);
    HistMean(i)= mean(LogRendFull,1);
    end
    else
    HistSigma = std(LogRend,1);
    HistMean= mean(LogRend,1);
    end
    CorrMat = corr(LogRend);
    CorrMat(NumAsset+1:end,:) = [];
    CorrMat(:,1:NumAsset) = [];
    HistRho = diag(CorrMat);
    HistRho(isnan(HistRho)) = 0;
    
    Pxt0 = zeros(1,NumAsset);
    for i=1:NumAsset
    Data1=Data(Data(:,i*2-1)>0,i*2);
    Pxt0(1,i) = Data1(end,1);
    end
    end