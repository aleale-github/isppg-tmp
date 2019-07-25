function [RiskPRIIPs,PerformancePRIIPs,CostPRIIPs,PxProdPerf,PxProdPerfStress,PxUndPerf,PxUndStress] = PRIIPsModelFramework(LogRend,Underlying,LogRendAsset,UnderlyingName,Freq,IRC,Div,FixingDate,PaymentFloatDate,PaymentFixDate,PricingDate,BasisCouponFreq,ParaProd,Inv,PxProd0,DistrFee,ManagFee,RunningFee,ExitFee,FixingType,ProductType,UnderlyingType,PathDep,Rating,Cat,PercPerf,NumSim,SimType,Seed,Figure)

    warning off
% From Excel to Matlab
    [IRC,FixingDate,PaymentFloatDate,PaymentFixDate,PricingDate] = DateMatWrkDays(IRC,FixingDate,PaymentFloatDate,PaymentFixDate,PricingDate);

% Date Handling
    [FixingDateFut,IssueDate,PaymentFloatDateFut,PaymentFixDateFut,PaymentDateInt,FixingUndPast0,PaymentDateIssue,PosDataNumHistSim] = DateHandling(Underlying,FixingDate,PaymentFloatDate,PaymentFixDate,PricingDate,[5;10],ProductType);

% Basic Underlying Analysis    
    [HistMean,Sigma,Rho,NumAsset,PxUnd0] = DataStat(LogRend,LogRendAsset,UnderlyingType); 
    
% Discount Factor into PRIIP currency
    Fwd = 0;
    [~,PtTFloat] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,end)],PricingDate,PaymentFloatDateFut,Fwd);
    [~,PtTFix] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,end)],PricingDate,PaymentFixDateFut,Fwd);
    [~,PtTInt] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,end)],PricingDate,PaymentDateInt,Fwd);

     if Cat == 2
    % Time Horizon
    %    PaymentDateFut = sort([PaymentFloatDateFut;PaymentFixDateFut]);    
    %    [N,~,T] = TimeSteps(PaymentDateFut(1),PaymentDateFut(2:end,1));               
    %    ParaCF = CornishFisherParaPRIIPs(LogRend(:,1:NumAsset));
    %    Scope = 'Risk';
    %    StressScen = 'N';
    %    PxProdRisk = CornishFisherPRIIPs(0.025,N(end),T(end),ParaCF,0,Scope,StressScen);
    %    vStress = StressVolatilityPRIIPs(LogRend(:,1:NumAsset),T(end),Freq);    
    %    Scope = 'Perf';
    %    StressScen = 'Y';
    %    PxProdPerf = CornishFisherPRIIPs([0.10;0.50;0.90],N,T,ParaCF,vStress,Scope,StressScen);
        
    % -------------------- RISK ANALYSIS -------------------- %
    %    [SRI,MRM,CRM,VEV] = RiskIndicatorPRIIPs(PxProdRisk,T(end),Rating,Cat);
        
    % ---------------- PERFORMANCE ANALYSIS ---------------- %        
    %    [RendNetFee,InvNetFee] = PerformanceScenariosPRIIPs(PxProdPerf,PxProd0,Inv,T,PaymentDateFut,UpFrontFee,RunningFee,ExitFee,PercPerf,Cat);
                
    elseif Cat == 3
    % Risk Parameters
        RiskFactor = find(cell2mat(UnderlyingType(2,:)) == 1);
        NumRiskFactor = size(RiskFactor,2);
        Fwd = cell2mat(UnderlyingType(end-2,RiskFactor));
        ParallelShift = cell2mat(UnderlyingType(end-1,RiskFactor));
    % Time Horizon
        [SimDateFut,PosDateFixing,PosDateInt,PosDatePayFloat,PaymentDate,PaymentDateFut,WrkDays,NumWrkDays,CouponFreq,T] = TimeData(PricingDate,PaymentFloatDateFut,PaymentFixDateFut,PaymentFloatDate,PaymentFixDate,PaymentDateInt,FixingDateFut,PathDep,BasisCouponFreq,ProductType);
    % Bootstrapping return with historical drift (Step 1)
        StressScen = 'N';    
        [Meanshift,RendUnd0Drift,RendUndHist]= ScenSimPRIIPs(LogRend(:,1:NumAsset),WrkDays,Sigma,StressScen,SimType,NumSim,Seed);
    % Bootstrapping with historical drift and stress volatility
        [SigmaStress,Stat] = StressVolatilityPRIIPs(LogRendAsset,LogRend(:,1:NumAsset),T(end),Freq(1,1),UnderlyingType);
        [SigmaStress1Y,Stat1Y] = StressVolatilityPRIIPs(LogRendAsset,LogRend(:,1:NumAsset),1,Freq(1,1),UnderlyingType);        
        StressScen = 'Y';
        [MeanshiftStress,RendUndStress0Drift,~] = ScenSimPRIIPs(LogRend(:,1:NumAsset),WrkDays,SigmaStress(1,1:NumAsset),StressScen,SimType,NumSim,Seed);
        [MeanshiftStress1Y,RendUndStress1Y0Drift,~]  = ScenSimPRIIPs(LogRend(:,1:NumAsset),WrkDays,SigmaStress1Y(1,1:NumAsset),StressScen,SimType,NumSim,Seed);

    % ---------------------------------------------------------- RISK CALCULATION ----------------------------------------------------------- %

    % Correction for Risk Neutral Drift    
        Drift = 'RF';
        DriftFwd = zeros(1,NumRiskFactor,size(SimDateFut,1));
        PxUndRisk = zeros(NumSim,NumRiskFactor,size(SimDateFut,1));
        FixingUndPast = zeros(NumRiskFactor,size(FixingUndPast0,1));        
        for j = 1:NumRiskFactor
            if strcmp(UnderlyingType{1,j},'Rate')
            % Fwd calculation    
                [~,~,rtTT] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,1+RiskFactor(j))],PricingDate,SimDateFut,Fwd(j));
                DriftFwd(1,1,:) = rtTT';
            % Step 3: correction of historical drift with Risk Neutral drift
                PxUndRisk1 = Ret2LevPRIIPs(Meanshift(:,RiskFactor(j)),RendUnd0Drift(:,RiskFactor(j),:),PxUnd0(1,RiskFactor(j)),ParallelShift);
                PxUndRisk1 = PxUndRisk1-repmat(mean(PxUndRisk1,1),NumSim,1,1)+repmat(DriftFwd,[NumSim,1,1]);
            elseif strcmp(UnderlyingType{1,j},'Asset')
            % Spot Calculation    
                [rtT] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,1+RiskFactor(j))],PricingDate,SimDateFut,Fwd(j));
            % Correction of historical drift with Risk Neutral drift
                RiskPremium = zeros(NumWrkDays,1);            
                rfree = rtT;
                d = repmat(Div(j)',NumWrkDays,1);
                sigmaUnd = repmat(Sigma(1,j),NumWrkDays,1);
                sigmaCurr = repmat(Sigma(1,j+RiskFactor(end)),NumWrkDays,1);
                Rho1=repmat(Rho(j,:),NumWrkDays,1); 

                if strcmp(PathDep{1},'D')
                    WrkDaysR = 365;
                elseif strcmp(PathDep{1},'W')
                    WrkDaysR = 52;
                elseif strcmp(PathDep{1},'M')
                    WrkDaysR = 12;
                end
                
                [RendUndRF] = DriftCorrection(RendUnd0Drift(:,RiskFactor(j),:),RiskPremium,rfree,d,sigmaUnd,sigmaCurr,Rho1,WrkDaysR,WrkDays,Drift);
                PxUndRisk1 = Ret2LevPRIIPs(Meanshift(:,RiskFactor(j)),RendUndRF,PxUnd0(1,RiskFactor(j)),ParallelShift(j));

            elseif strcmp(UnderlyingType{1,j},'Currency')
            % Spot Calculation    
                [rtTDiv1] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,1+RiskFactor(j))],PricingDate,SimDateFut,Fwd(j));
                [rtTDiv2] = DiscountFactorPRIIPs([IRC(:,1),IRC(:,end)],PricingDate,SimDateFut,Fwd(j));
            % Correction of historical drift with Risk Neutral drift
                RiskPremium = zeros(NumWrkDays,1);            
                rfree = rtTDiv1;
                d = rtTDiv2;
                sigmaUnd = repmat(Sigma(1,j),NumWrkDays,1);
                sigmaCurr = repmat(Sigma(1,j+RiskFactor(end)),NumWrkDays,1);
                Rho1=repmat(Rho(j,:),NumWrkDays,1); 
                if strcmp(PathDep{1},'D')
                    WrkDaysR = 365;
                elseif strcmp(PathDep{1},'W')
                    WrkDaysR = 52;
                elseif strcmp(PathDep{1},'M')
                    WrkDaysR = 12;
                end
                [RendUndRF] = DriftCorrection(RendUnd0Drift(:,RiskFactor(j),:),RiskPremium,rfree,d,sigmaUnd,sigmaCurr,Rho1,WrkDaysR,WrkDays,Drift);
                PxUndRisk1 = Ret2LevPRIIPs(Meanshift(:,RiskFactor(j)),RendUndRF,PxUnd0(1,RiskFactor(j)),ParallelShift(j));
            end
            PxUndRisk(:,j,:) = PxUndRisk1;
            FixingUndPast(j,:) = FixingUndPast0(:,RiskFactor(j));
        end
            FixingUndPast = FixingUndPast';
            
    % Underlying Px Scenarios
        if strcmp(ProductType,'TC') || strcmp(ProductType,'TCA') || ... 
            strcmp(ProductType,'TCS') || strcmp(ProductType,'TCSA') || ...
            strcmp(ProductType,'TPPF') || strcmp(ProductType,'TPPFA') || ...
            strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPUA') || ... 
            strcmp(ProductType,'TPPVI') || ...
            strcmp(ProductType,'CPPVI') || strcmp(ProductType,'CPPVE') || ...
            strcmp(ProductType,'CCTI') || strcmp(ProductType,'CCTE') || ...
            strcmp(ProductType,'CPBFI') || strcmp(ProductType,'CPBFE') || ...
            strcmp(ProductType,'CPFKII') || strcmp(ProductType,'CPFKIE') || ...
            strcmp(ProductType,'CPFKIKOI') || strcmp(ProductType,'CPFKIKOE') || ...
            strcmp(ProductType,'CPCI') || strcmp(ProductType,'CPCE') || ...
            strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV') || ...
            strcmp(ProductType,'CPCA2') || strcmp(ProductType,'CPCV2') || ...
            strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE') || ...
            strcmp(ProductType,'CDAI') || strcmp(ProductType,'CDAE') || ...
            strcmp(ProductType,'CDALI') || strcmp(ProductType,'CDALE')
            Method = 'DRM2';
        else
            Method = 'DRM3';
        end
        
        %{
        if strcmp(ProductType,'FRC')
            PxProd0Val = PxProd0;
            DistrFeeVal = DistrFee;
            ManagFeeVal = ManagFee;
            RunningFeeVal = RunningFee;
            ExitFeeVal = ExitFee;    
            if PricingDate >= IssueDate
                FlagCCSSim = 0;
                PxProd0 = PxProd0Val./Underlying(Underlying(:,1)==PricingDate,2);
                DistrFee = DistrFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ManagFee = ManagFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                RunningFee = RunningFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ExitFee = ExitFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
            elseif PricingDate < IssueDate
                FlagCCSSim = 1;
                PxProd0 = PxProd0Val./PxUndRisk(:,1,1);
                DistrFee = DistrFeeVal./PxUndRisk(:,1,1);
                ManagFee = ManagFeeVal./PxUndRisk(:,1,1);
                RunningFee = RunningFeeVal./PxUndRisk(:,1,1);
                ExitFee = ExitFeeVal./PxUndRisk(:,1,1);
            end
        end
        %}
        
        Scope = 'Risk'; Kernel = 0;
        [PxProdRisk,~,~,PriceDC] = ProductPricing(PxProd0,ParaProd,PxUndRisk,PtTFloat,PtTFix,PtTInt,PosDateFixing,PosDatePayFloat,PosDateInt,FixingUndPast,PricingDate,PaymentDateFut,PaymentDateInt,PaymentDateIssue,IssueDate,Method,CouponFreq,FixingType,UnderlyingType,ProductType,PathDep,Scope,Kernel);
        if strcmp(ProductType,'FRC')
         PxProdRisk=PxProdRisk.*PriceDC;
         end
        [SRI,MRM,CRM,VEV] = RiskIndicatorPRIIPs(PxProd0,PxProdRisk,T(end),DistrFee,Rating,Cat);
        
    % ------------------------------------------------------- PERFORMANCE CALCULATION ------------------------------------------------------- %

    % Correction for Historical Drift
        Drift = 'HR';
        PxUndPerf = zeros(NumSim,NumRiskFactor,size(SimDateFut,1));
        PxUndStress = zeros(NumSim,NumRiskFactor,size(SimDateFut,1));
        PxUndStress1Y = zeros(NumSim,NumRiskFactor,size(SimDateFut,1));
    % Setting Parameters
        for j = 1:NumRiskFactor
            RiskPremium = repmat(HistMean(1,RiskFactor(j)),NumWrkDays,1);
            rfree = zeros(NumWrkDays,RiskFactor(j));
            d = zeros(NumWrkDays,RiskFactor(j));
            sigmaUnd = repmat(Sigma(1,RiskFactor(j)),NumWrkDays,1);
            sigmaCurr = zeros(NumWrkDays,1);
            if strcmp(PathDep{1},'D')
                WrkDaysR = 365;
            elseif strcmp(PathDep{1},'W')
                WrkDaysR = 52;
            elseif strcmp(PathDep{1},'M')
                WrkDaysR = 12;
            end
        % Performance Base Prices
            %[RendUndHist] = DriftCorrection(RendUnd0Drift(:,RiskFactor(j),:),RiskPremium,rfree,d,sigmaUnd,sigmaCurr,Rho,WrkDaysR,WrkDays,Drift);
            PxUndPerf1 = Ret2LevPRIIPs(Meanshift(:,RiskFactor(j)),RendUndHist(:,RiskFactor(j),:),PxUnd0(1,RiskFactor(j)),ParallelShift(j));
            %PxUndPerf1 = Ret2LevPRIIPs(RendUndHist(:,1,:),PxUnd0(1,RiskFactor(j)),ParallelShift(j));
            PxUndPerf(:,j,:) = PxUndPerf1;
        % Performance Stress > 1Y Prices
            sigmaUndStress = repmat(SigmaStress(1,RiskFactor(j)),NumWrkDays,1);            
            RiskPremium = -0.5*sigmaUndStress.^2;
            [RendUndStress0DriftVol] = DriftCorrection(RendUndStress0Drift(:,RiskFactor(j),:),RiskPremium,rfree,d,sigmaUndStress,sigmaCurr,Rho,WrkDaysR,WrkDays,Drift);
            PxUndStress1 = Ret2LevPRIIPs(MeanshiftStress(:,RiskFactor(j)),RendUndStress0DriftVol,PxUnd0(1,RiskFactor(j)),ParallelShift(j));
            PxUndStress(:,j,:) = PxUndStress1;
        % Performance Stress <= 1Y Prices
            sigmaUndStress1Y = repmat(SigmaStress1Y(1,RiskFactor(j)),NumWrkDays,1);            
            RiskPremium1Y = -0.5*sigmaUndStress1Y.^2;            
            [RendUndStress1Y0DriftVol] = DriftCorrection(RendUndStress1Y0Drift(:,RiskFactor(j),:),RiskPremium1Y,rfree,d,sigmaUndStress1Y,sigmaCurr,Rho,WrkDaysR,WrkDays,Drift);
            PxUndStress1Y1 = Ret2LevPRIIPs(MeanshiftStress1Y(:,RiskFactor(j)),RendUndStress1Y0DriftVol,PxUnd0(1,RiskFactor(j)),ParallelShift(j));
            PxUndStress1Y(:,j,:) = PxUndStress1Y1;
        end            
    % Underlying Px Scenarios
        PricingDate = max(PricingDate,PaymentFixDate(1));
        Scope = 'Perf';
        %{
        if strcmp(ProductType,'FRC')
            PxProd0Val = PxProd0;
            DistrFeeVal = DistrFee;
            ManagFeeVal = ManagFee;
            RunningFeeVal = RunningFee;
            ExitFeeVal = ExitFee;    
            if FlagCCSSim == 0
                PxProd0 = PxProd0Val./Underlying(Underlying(:,1)==PricingDate,2);
                DistrFee = DistrFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ManagFee = ManagFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                RunningFee = RunningFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ExitFee = ExitFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
            elseif FlagCCSSim == 1
                PxProd0 = PxProd0Val./PxUndPerf(:,1,1);
                DistrFee = DistrFeeVal./PxUndPerf(:,1,1);
                ManagFee = ManagFeeVal./PxUndPerf(:,1,1);
                RunningFee = RunningFeeVal./PxUndPerf(:,1,1);
                ExitFee = ExitFeeVal./PxUndPerf(:,1,1);
            end
        end
        %}
        [~,PxProdPerf,PxProdPerfCF,PriceDC,DFCostPerf] = ProductPricing(PxProd0+DistrFee,ParaProd,PxUndPerf,PtTFloat,PtTFix,PtTInt,PosDateFixing,PosDatePayFloat,PosDateInt,FixingUndPast,PricingDate,PaymentDateFut,PaymentDateInt,PaymentDateIssue,IssueDate,Method,CouponFreq,FixingType,UnderlyingType(:,RiskFactor),ProductType,PathDep,Scope,Kernel);
    % Performance Scenarios (Normal)
        [RendNetFee3,InvNetFee3,RendCost,InvCost] = PerformanceCostScenariosPRIIPs(PriceDC,PxProdPerf,PxProdPerfCF,DFCostPerf,PxProd0,ParaProd,Inv,T,FixingDate,PaymentDate,PaymentFixDateFut,PaymentDateInt,PosDateInt,PosDateFixing,DistrFee,ManagFee,RunningFee,ExitFee,PercPerf,Cat,Method,ProductType,Scope);
    % Performance Scenarios
        Scope = 'Stress';
    % Stress > 1Y
    %{
        if strcmp(ProductType,'FRC')
            PxProd0Val = PxProd0;
            DistrFeeVal = DistrFee;
            ManagFeeVal = ManagFee;
            RunningFeeVal = RunningFee;
            ExitFeeVal = ExitFee;    
            if FlagCCSSim == 0
                PxProd0 = PxProd0Val./Underlying(Underlying(:,1)==PricingDate,2);
                DistrFee = DistrFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ManagFee = ManagFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                RunningFee = RunningFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ExitFee = ExitFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
            elseif FlagCCSSim == 1
                PxProd0 = PxProd0Val./PxUndStress(:,1,1);
                DistrFee = DistrFeeVal./PxUndStress(:,1,1);
                ManagFee = ManagFeeVal./PxUndStress(:,1,1);
                RunningFee = RunningFeeVal./PxUndStress(:,1,1);
                ExitFee = ExitFeeVal./PxUndStress(:,1,1);
            end
        end
    %}
        [~,PxProdPerfStress,PxProdPerfStressCF,PriceDCStress,DFCostPerf] = ProductPricing(PxProd0+DistrFee,ParaProd,PxUndStress,PtTFloat,PtTFix,PtTInt,PosDateFixing,PosDatePayFloat,PosDateInt,FixingUndPast,PricingDate,PaymentDateFut,PaymentDateInt,PaymentDateIssue,IssueDate,Method,CouponFreq,FixingType,UnderlyingType(:,RiskFactor),ProductType,PathDep,Scope,Kernel);
        [RendNetFeeStress,InvNetFeeStress] = PerformanceCostScenariosPRIIPs(PriceDCStress,PxProdPerfStress,PxProdPerfStressCF,DFCostPerf,PxProd0,ParaProd,Inv,T,FixingDate,PaymentDate,PaymentFixDateFut,PaymentDateInt,PosDateInt,PosDateFixing,DistrFee,ManagFee,RunningFee,ExitFee,5,Cat,Method,ProductType,Scope);
    % Stress <= 1Y
    %{
        if strcmp(ProductType,'FRC')
            PxProd0Val = PxProd0;
            DistrFeeVal = DistrFee;
            ManagFeeVal = ManagFee;
            RunningFeeVal = RunningFee;
            ExitFeeVal = ExitFee;    
            if FlagCCSSim == 0
                PxProd0 = PxProd0Val./Underlying(Underlying(:,1)==PricingDate,2);
                DistrFee = DistrFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ManagFee = ManagFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                RunningFee = RunningFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
                ExitFee = ExitFeeVal./Underlying(Underlying(:,1)==PricingDate,2);
            elseif FlagCCSSim == 1
                PxProd0 = PxProd0Val./PxUndStress1Y(:,1,1);
                DistrFee = DistrFeeVal./PxUndStress1Y(:,1,1);
                ManagFee = ManagFeeVal./PxUndStress1Y(:,1,1);
                RunningFee = RunningFeeVal./PxUndStress1Y(:,1,1);
                ExitFee = ExitFeeVal./PxUndStress1Y(:,1,1);
            end
        end
    %}
        [~,PxProdPerfStress1Y,PxProdPerfStressCF1Y,PriceDCStress1Y,DFCostPerf1Y] = ProductPricing(PxProd0+DistrFee,ParaProd,PxUndStress1Y,PtTFloat,PtTFix,PtTInt,PosDateFixing,PosDatePayFloat,PosDateInt,FixingUndPast,PricingDate,PaymentDateFut,PaymentDateInt,PaymentDateIssue,IssueDate,Method,CouponFreq,FixingType,UnderlyingType(:,RiskFactor),ProductType,PathDep,Scope,Kernel);
        [RendNetFeeStress1Y,InvNetFeeStress1Y] = PerformanceCostScenariosPRIIPs(PriceDCStress1Y,PxProdPerfStress1Y,PxProdPerfStressCF1Y,DFCostPerf1Y,PxProd0,ParaProd,Inv,T,FixingDate,PaymentDate,PaymentFixDateFut,PaymentDateInt,PosDateInt,PosDateFixing,DistrFee,ManagFee,RunningFee,ExitFee,5,Cat,Method,ProductType,Scope);

        RendNetFee = [RendNetFeeStress;RendNetFee3];
        InvNetFee = [InvNetFeeStress;InvNetFee3];
        if T(end) > 1
            RendNetFee(1,1) = RendNetFeeStress1Y(1);
            InvNetFee(1,1) = InvNetFeeStress1Y(1);
        end
    end
    
% OUTPUT MANAGEMENT
    NumPerfScen = size(RendNetFee,1);
    NumHP = size(RendNetFee,2);
    PerformancePRIIPs = zeros(NumPerfScen*2,NumHP);
    CostPRIIPs = [InvCost;RendCost];
    p = (1:2:2*NumPerfScen);
    for i = 1:NumPerfScen
        PerformancePRIIPs(p(i),:) = InvNetFee(i,:);
        PerformancePRIIPs(p(i)+1,:) = RendNetFee(i,:);
    end
    RiskPRIIPs = [VEV;MRM;CRM;SRI];

% Additional Statistics    

% FIGURE
    if strcmp(Figure{1},'Y')
    % historical Series
        figure
        for i = 1:NumRiskFactor    
            subplot(NumRiskFactor,1,i)
            Underlying1=Underlying(Underlying(:,RiskFactor(i)*2-1)>0,RiskFactor(i)*2-1:RiskFactor(i)*2);
            plot(Underlying1(end-PosDataNumHistSim(i)+1:end,1),Underlying1(end-PosDataNumHistSim(i)+1:end,2));
            a = gca;
            a.XTick = linspace(Underlying1(end-PosDataNumHistSim(1)+1,1),Underlying1(end,1),6);
            grid on
            title([UnderlyingName(RiskFactor(i)),' Historical Series'])
            datetick('x','mmmyy','keepticks') 
            xlabel('Date')
            ylabel('Level')
            if strcmp(UnderlyingType(1,RiskFactor(i)),'Rate') 
                a.YAxis.TickLabelFormat = '%.2F%%';                
            end
        end
        saveas(gcf,'C:\Users\carlomaria.petrangel\Desktop\PRIIPS\Figures\Figure1','bmp');

    % Historical Log Distribution
        figure
        for i = 1:NumRiskFactor
            subplot(NumRiskFactor,1,i)
            Underlying1=Underlying(Underlying(:,RiskFactor(i)*2-1)>0,RiskFactor(i)*2-1:RiskFactor(i)*2);
            plot(Underlying1(end-size(LogRend,1)+1:end,1),LogRend(:,RiskFactor(i)).*100);
            a = gca;
            a.XTick = linspace(Underlying1(end-size(LogRend,1)+1,1),Underlying1(end,1),6);            
            grid on
            title([UnderlyingName(RiskFactor(i)),' Log Return'])
            datetick('x','mmmyy','keepticks') 
            xlabel('Date')
            ylabel('Log Return')
            a.YAxis.TickLabelFormat = '%.2F%%';
            legend(UnderlyingName(RiskFactor(i)),'location','best')
        end            
        saveas(gcf,'C:\Users\carlomaria.petrangel\Desktop\PRIIPS\Figures\Figure2','bmp');

    % Historical Underlying Simulation
        figure
        NumScen = size(PxProdPerf,2);        
        for i = 1:NumRiskFactor
            subplot(NumRiskFactor,1,i)
            t = round(WrkDays./252,4);
            plot(t,prctile(squeeze(PxUndPerf(:,i,:)),[1 10 20 50 80 90 99])','o-');
            a = gca;
            if size(T,1) == 1
                a.XTick = round(T,2);
            else
                %a.XTick = round(linspace(T(1),T(end,1),NumScen*2),2);
                a.XTick = round(linspace(T(1),T(end,1),NumScen),2);
            end            
            grid on
            title([UnderlyingName(RiskFactor(i)),' Px Simulation'])
            xlabel('Simulation Steps (Days)')
            ylabel('Und Level')
            legend('Perc_{1%}','Perc_{10%}','Perc_{20%}','Perc_{50%}','Perc_{80%}','Perc_{90%}','Perc_{99%}','location','best')
            if strcmp(UnderlyingType(1,RiskFactor(i)),'Rate') 
                a.YAxis.TickLabelFormat = '%.2F%%';                
            end
        end
        saveas(gcf,'C:\Users\carlomaria.petrangel\Desktop\PRIIPS\Figures\Figure3','bmp');
        
    % Historical Stress Volatility
        figure
        for i = 1:NumRiskFactor
            subplot(NumRiskFactor,1,i)
            NumVolStress1Y = size(Stat1Y,1);
            Underlying1=Underlying(Underlying(:,RiskFactor(i)*2-1)>0,RiskFactor(i)*2-1:RiskFactor(i)*2);
            if T(end) > 1
                NumVolStress = size(Stat,1);
                Stat1 = Stat(:,RiskFactor(i));
                plot(Underlying1(end-NumVolStress+1:end,1),Stat1.*sqrt(252).*100);
                hold on
                Stat1Y1 = [Stat1Y(:,RiskFactor(i)),repmat(Sigma(1,RiskFactor(i)),NumVolStress1Y,1),repmat(SigmaStress(1,RiskFactor(i)),NumVolStress1Y,1),repmat(SigmaStress1Y(1,RiskFactor(i)),NumVolStress1Y,1)];
            else
                Stat1Y1 = [Stat1Y(:,RiskFactor(i)),repmat(Sigma(1,RiskFactor(i)),NumVolStress1Y,1),repmat(SigmaStress1Y(1,RiskFactor(i)),NumVolStress1Y,1)];
            end
            plot(Underlying1(end-NumVolStress1Y+1:end,1),Stat1Y1.*sqrt(252).*100);
            hold off
            a = gca;
            a.XTick = linspace(Underlying1(end-NumVolStress1Y+1,1),Underlying1(end,1),6);            
            grid on
            title('Historical Stress Volatility')
            datetick('x','mmmyy','keepticks') 
            xlabel('Date')
            ylabel('$$^{w}_{t_{i}}\sigma_{Stress} = \frac{\sum_{t_{i}}^{t_{i+w}}{(r_{i}-^{t{i+w}}_{t_{i}}M_{1})^{2}}} {M_{w}}$$','Interpreter','latex')
            a.YAxis.TickLabelFormat = '%.2F%%';
            if T(end) > 1
                legend('\sigma^{Hist}_{Stress > 1Y}','\sigma^{Hist}_{Stress \leq 1Y}','\sigma^{Perc}_{Base}','\sigma^{Perc}_{Stress > 1Y}','\sigma^{Perc}_{Stress \leq 1Y}','location','best')
            else
                legend('\sigma^{Hist}_{Stress \leq 1Y}','\sigma^{Perc}_{Base}','\sigma^{Perc}_{Stress \leq 1Y}','location','best')
            end
        end
        saveas(gcf,'C:\Users\carlomaria.petrangel\Desktop\PRIIPS\Figures\Figure5','bmp');
        
     % Historical Product Distribution
        figure
        NumScen = size(PxProdPerf,2);
        for t = 1:NumScen
        if strcmp(ProductType,'TC') || strcmp(ProductType,'TCA') || ...
            strcmp(ProductType,'TCS') || strcmp(ProductType,'TCSA') || ...
            strcmp(ProductType,'TPPF') || strcmp(ProductType,'TPPFA') || ...
            strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPUA') || ... 
            strcmp(ProductType,'TPPVI') || ...
            strcmp(ProductType,'CPPVI') || strcmp(ProductType,'CPPVE') || ...
            strcmp(ProductType,'CCTI') || strcmp(ProductType,'CCTE') || ...
            strcmp(ProductType,'CPFKII') || strcmp(ProductType,'CPFKIE') || ...
            strcmp(ProductType,'CPFKIKOI') || strcmp(ProductType,'CPFKIKOE') || ...
            strcmp(ProductType,'CPCI') || strcmp(ProductType,'CPCE') || ...
            strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV') || ...
            strcmp(ProductType,'CPCA2') || strcmp(ProductType,'CPCV2') || ...
            strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE') || ...
            strcmp(ProductType,'CDAI') || strcmp(ProductType,'CDAE') || ...
            strcmp(ProductType,'CDALI') || strcmp(ProductType,'CDALE')
            PDF = PxProdPerf{:,t};
        else
            PDF = PxProdPerf{:,t};
        end
            plot(T(t),prctile(PDF,[1 10 20 50 80 90 99]) ,'^--','MarkerFaceColor','auto','MarkerSize',10);
            hold on
        end
        a = gca;
        if size(T,1) == 1
            a.XTick = round(T,2);
        else
            %a.XTick = round(linspace(T(1),T(end,1),NumScen*2),2);   
            a.XTick = round(linspace(T(1),T(end,1),NumScen),2);
        end
        grid on
        title('Product Px Simulation')
        xlabel('Simulation Steps (Days)')
        ylabel('Product Price')
        legend('Perc_{1%}','Perc_{10%}','Perc_{20%}','Perc_{50%}','Perc_{80%}','Perc_{90%}','Perc_{99%}','location','northwest')
        if strcmp(UnderlyingType{1,1},'Px') 
            a.YAxis.TickLabelFormat = '%.2F%%';                
        end
        saveas(gcf,'C:\Users\carlomaria.petrangel\Desktop\PRIIPS\Figures\Figure4','bmp');             
        
    % PDF Perf Base vs PDF Perf Stress
        figure
        a = gca;
        if strcmp(ProductType,'TC') || strcmp(ProductType,'TCA') || ...
            strcmp(ProductType,'TCS') || strcmp(ProductType,'TCSA') || ...
            strcmp(ProductType,'TPPF') || strcmp(ProductType,'TPPFA') || ...
            strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPUA') || ... 
            strcmp(ProductType,'TPPVI') || ...
            strcmp(ProductType,'CPPVI') || strcmp(ProductType,'CPPVE') || ...
            strcmp(ProductType,'CCTI') || strcmp(ProductType,'CCTE') || ...
            strcmp(ProductType,'CPFKII') || strcmp(ProductType,'CPFKIE') || ...
            strcmp(ProductType,'CPFKIKOI') || strcmp(ProductType,'CPFKIKOE') || ...
            strcmp(ProductType,'CPCI') || strcmp(ProductType,'CPCE') || ...
            strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')|| ...
            strcmp(ProductType,'CPCA2') || strcmp(ProductType,'CPCV2') || ...
            strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE') || ...
            strcmp(ProductType,'CDAI') || strcmp(ProductType,'CDAE') || ...
            strcmp(ProductType,'CDALI') || strcmp(ProductType,'CDALE')
            histogram(PxProdPerf{:,end},50)
            hold on
            histogram(PxProdPerfStress{:,end},50)
            if T(end) > 1
                histogram(PxProdPerfStress1Y{:,1},50)
            end
            grid on
        else        
            histogram(PxProdPerf{:,end}*100,50)
            hold on
            histogram(PxProdPerfStress{:,end}*100,50)
            if T(end) > 1
                histogram(PxProdPerfStress1Y{:,1}*100,50)
            end
            grid on
        end
        title('PDF Product Price Performance on RHP: Base vs Stress')
        xlabel('Price')
        a.XAxis.TickLabelFormat = '€%.2f';
        ylabel('F(x)')
        if T(end) > 1
            legend('PDF^{Perf}_{Base}','PDF^{Perf}_{Stress}','PDF^{Perf}_{Stress1Y}','location','northwest')
        else
            legend('PDF^{Perf}_{Base}','PDF^{Perf}_{Stress}','location','northwest')
        end
        saveas(gcf,'C:\Users\carlomaria.petrangel\Desktop\PRIIPS\Figures\Figure6','bmp');
   close all
    end
