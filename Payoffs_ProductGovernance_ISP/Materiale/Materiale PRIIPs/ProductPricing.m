function [PxProdRisk,PxProdPerf,PxProdPerfCF,PriceDC,DFCostPerf] = ProductPricing(PxProd0,ParaProd,PxUndFloat,PtTFloat,PtTFix,PtTInt,PosDateFixing,PosDatePayFloat,PosDateInt,FixingUndPast,PricingDate,PaymentDateFut,PaymentDateInt,PaymentDateIssue,IssueDate,Method,CouponFreq,FixingType,DataType,ProductType,PathDep,Scope,Kernel)
% Pricing Function    

% Internal Variable
    NumScenInt = size(PtTInt,1);
    NumSim = size(PxUndFloat,1);
    PxProdPerf = cell(1,NumScenInt+1);
    PxProdPerfCF = cell(1,NumScenInt+1);    
    NumCouponFixIssue = size(PaymentDateIssue,1);

% Repricing on RHP
    if PaymentDateFut(1) ==  IssueDate
        PaymentDateFut(1) = [];
        PtTFloat(1,:) = [];
        PtTFix(1,:) = [];
        PosDatePayFloat = PosDatePayFloat-1;
        PosDatePayFloat(PosDatePayFloat == 0) = [];
    end
    NumPayDate = size(PtTFix,1);        
    if strcmp(Scope,'Risk')
        [PxProdT,~,CouponFix,PriceDC] = ProductPayoff(ParaProd,PxUndFloat,CouponFreq,PosDateFixing,PosDateInt,NumCouponFixIssue,FixingUndPast,FixingType,PathDep,ProductType,DataType);
    % Selection of scenarios > Pricing Date    
        PxProdT = PxProdT(:,end-size(PtTFloat,1)+1:end);
        CouponFix = CouponFix(:,end-NumPayDate+1:end);
    % Discounting Price Distribution in t0
        PxProdRisk = PxProdT*PtTFloat+CouponFix*PtTFix;

    elseif strcmp(Scope,'Perf') || strcmp(Scope,'Stress')
        PxProdT1 = zeros(NumSim,NumPayDate);        
        [PxProdT,PxUndFloat,CouponFix,PriceDC] = ProductPayoff(ParaProd,PxUndFloat,CouponFreq,PosDateFixing,PosDateInt,NumCouponFixIssue,FixingUndPast,FixingType,PathDep,ProductType,DataType);
    % Selection of scenarios > Pricing Date    
        PxProdT = PxProdT(:,end-size(PtTFloat,1)+1:end);
            CouponFix = CouponFix(:,end-NumPayDate+1:end);
    % Resizing Distribution on all Payment Date  
        PxProdT1(:,PosDatePayFloat) = PxProdT;
    % Capitalization Price Distribution in RHP         
        PxProdPerf{1,end} = ((PxProdT*PtTFloat)+(CouponFix*PtTFix))./PtTFix(end);
        
        if strcmp(ProductType,'FRC')
        PxProdPerf{1,end}=PxProdPerf{1,end}.*PriceDC;
        end
        
        DFCostPerf = cell(1,NumScenInt+1);
        DFCostPerf{1,end} = 1./PtTFix;
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
            strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')|| ...
            strcmp(ProductType,'CPCA2') || strcmp(ProductType,'CPCV2') || ...
            strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE') || ...
            strcmp(ProductType,'CDAI') || strcmp(ProductType,'CDAE') || ...
            strcmp(ProductType,'CDALI') || strcmp(ProductType,'CDALE')
            %PxProdPerf{:,end} = PxProdPerf{:,end}+1;
            PxProdPerf{:,end} = PxProdPerf{:,end}-PxProd0(:,1)./PtTFix(end);
            PxProdPerf{:,end} = PxProdPerf{:,end};            
            %DFCostPerf(1,end) = 1./PtTFix(end);
            PxProd0CF=zeros(size(PxProdT1,1),size(PxProdT1,2));
            PxProd0CF(:,1)=-PxProd0(:,1)./PtTFix(1);
            PxProdPerfCF{1,end} = PxProdT1+CouponFix+PxProd0CF;
            PxProdPerfCF1 = PxProdT1+CouponFix+PxProd0CF;
        else        
            PxProdPerfCF{1,end} = PxProdT1+CouponFix;
            PxProdPerfCF1 = PxProdT1+CouponFix;
        end
        % Repricing on Intermediate Date   
            X = @(x) [ones(size(x,1),1) x x.^2 x.^3 x.^4 x.^5];     % Function F(x): y = a0+a1*X+a2*X^2+a3*X^3+a4*X^4+a5*X^5
            %X = @(x) [exp(-x./2) -2.*x.*exp(-x./2) 0.5.*x.^2.*exp(-x./2)];
        if isempty(NumScenInt) == 0
            for IS = 1:NumScenInt
            % Underlying (X)    
                if strcmp(PathDep{1},'D') || strcmp(PathDep{1},'W') || strcmp(PathDep{1},'M')
                    PosDatePast = size(PaymentDateFut(PaymentDateFut(:,1) <= PaymentDateInt(IS),1),1);
                    if PosDatePast == 0
                        PosDatePast = 1;
                        PxProdPerfCFPast = [];
                        PxProdPerfCFFut = PxProdPerfCF1(:,PosDatePast:end);                
                    else
                        PxProdPerfCFPast = PxProdPerfCF1(:,1:PosDatePast);
                        PxProdPerfCFFut = PxProdPerfCF1(:,PosDatePast+1:end);
                    end
                    PxUndInt = PxUndFloat(:,:,IS);
                    NumCpnPastInt = size(PxProdPerfCFPast,2);
                    NumCpnFutInt = size(PxProdPerfCFFut,2);
                end
                s = squeeze(reshape(PxUndInt,size(PxUndInt,1),size(PxUndInt,2)*size(PxUndInt,3)));
                        
                if strcmp(ProductType,'FRC')
                    
                    PxProdPerfCFFut=PxProdPerfCFFut.*PriceDC;
                    PxProdPerfCFPast=PxProdPerfCFPast.*PriceDC;
                end
                
            % Payoff (Y)    
                if strcmp(Method,'DRM1') 
                    CashFlow = [repmat(-PxProd0,NumSim,1),PxProdPerf{:,end}];
                    CashFlow = mean(CashFlow,1)';                    
                    Date = [PricingDate;PaymentDateFut(end,1)];
                    DateInt = [PricingDate;PaymentDateInt];
                    dt = yearfrac(Date(1),Date(2:end),12);
                    dtInt = yearfrac(Date(1),DateInt(2:end),12);
                % Forward Risk Free floor calculation
                    PtTT = PtTFix(end)./PtTInt(IS);
                    rFree = (1/PtTT).^(1./(dt-dtInt(IS)))-1;
                    r = xirr(CashFlow,Date,0.01,1000,3);
                    r = max(rFree,r);
                    P0T = 1./((1+r).^dt(end));
                    P0T2 = 1./((1+r).^dtInt(IS));
                    PtTT = P0T./P0T2;
                    y = PxProdPerf{1,end}*PtTT;
                    PxProdPerfCF1 = PxProdPerfCF{1,end};
                elseif strcmp(Method,'DRM2')
                    PtTT = PtTFix./PtTInt(IS);
                    %y = (PxProdPerf{1,end})*PtTT;
                    y = PxProdPerfCFFut*PtTT(end-NumCpnFutInt+1:end,1);     
                    DFCostPerf1 = DFCostPerf{1,end};                        
                    %DFCostPerf{1,IS} = DFCostPerf1(1,end)*PtTT;
                    DFCostPerf{1,IS} = DFCostPerf1(1,end)*PtTT(end);
                    PxProdPerfCF1 = PxProdPerfCF{1,end};
                elseif strcmp(Method,'DRM3')
                    CashFlow = [repmat(-PxProd0,NumSim,1),PxProdPerf{:,end}];
                    CashFlow = mean(CashFlow,1)';                    
                    Date = [PricingDate;PaymentDateFut];
                    Date1 = [PricingDate;PaymentDateFut(end)];
                    DateInt = [PricingDate;PaymentDateInt];
                    dt = yearfrac(Date(1),Date(2:end),12);
                    dtInt = yearfrac(Date(1),DateInt(2:end),12); 
                % Forward Risk Free floor calculation
                    PtTT = PtTFix./PtTInt(IS);
                    rFree = (1./PtTT).^(1./(dt-dtInt(IS)))-1;
                    r = xirr(CashFlow,Date1,0.01,1000,3);
                    rStar = max(rFree(end),r);
                    %P0T = 1./((1+rStar).^dt(end-NumCpnFutInt+1:end,1));
                    %P0T2 = 1./((1+rStar).^dtInt(IS));
                    %PtTTStar = P0T./P0T2;
                    PtTTStar = 1./((1+rStar).^(dt-dtInt(IS)));
                    PtTT(end-NumCpnFutInt+1:end,1) = PtTTStar(end-NumCpnFutInt+1:end,1);
                    y = PxProdPerfCFFut*PtTT(end-NumCpnFutInt+1:end,1);                    
                end
                               
                if IS >= 1
                    if strcmp(ProductType,'ATC1') || strcmp(ProductType,'ATC2') || strcmp(ProductType,'ATD1') || strcmp(ProductType,'ATD2') || strcmp(ProductType,'ATD3') || strcmp(ProductType,'WOATD') ||strcmp(ProductType,'WOATD1')|| strcmp(ProductType,'ATCTW') 
                        if PosDatePast == 0
                            Aut = zeros(NumSim,1);
                        else
                            Aut = sum(PxProdPerfCF1(:,1:PosDatePast),2) >= 1;
                        end
                    else
                        Aut = zeros(NumSim,1);
                    end
                end
                PxProdt = y;
                s(Aut == 1,:) = [];
                y(Aut == 1,:) = [];
                if Kernel == 0
                    if PosDateInt(IS) >= PosDateFixing(end)
                        %PxProdt1 = y./PtTT./PtTFix(end).*PtTInt(IS);
                        PxProdt1 = y;
                    else
                        Coef = (X(s)'*X(s))\X(s)'*y;
                        PxProdt1 = (Coef'*X(s)')';
                    end
                elseif Kernel == 1
                    delta = 6;
                    PxProdt1 = NonParRegressionNorm(s',y',s',delta)';
                end
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
                    PxProdt(Aut == 0) = PxProdt1;                                           
                else
                    PxProdt(Aut == 0) = PxProdt1;
                end               
                if strcmp(Method,'DRM1')          
                    PxProdPerf{:,IS} = PxProdt;
                    PxProdPerfCFPast = PxProdPerfCF{1,end};
                    if PosDatePast == 0
                        PxProdPerfCFPast = zeros(NumSim,1);
                    else
                        PxProdPerfCFPast = PxProdPerfCFPast(:,1:PosDatePast);
                    end
                    PxProdPerfCF{1,IS} = PxProdPerfCFPast;                    
                elseif strcmp(Method,'DRM3')
                    if NumCpnPastInt == 0
                        PxProdPerf{:,IS} = PxProdt;
                    else
                        % PxProdPerf{:,IS} = PxProdPerfCFPast*(1./PtTT(1:NumCpnPastInt,1))+PxProdt;
                        PxProdPerf{:,IS} = PxProdPerfCFPast*(PtTT(1:NumCpnPastInt,1))+PxProdt;
                    end
                    PxProdPerfCF{1,IS} = PxProdPerfCFPast;
                elseif strcmp(Method,'DRM2')
                  if NumCpnPastInt == 0
                        PxProdPerf{:,IS} = PxProdt;
                  else
                        PxProdPerf{:,IS} = PxProdPerfCFPast*(PtTT(1:NumCpnPastInt,1))+PxProdt;
                        % PxProdPerf{:,IS} = PxProdPerfCFPast*(1./PtTT(1:NumCpnPastInt,1))+PxProdt;
                  end
                    PxProdPerfCF{1,IS} = PxProdPerfCFPast;  
%                     PxProdPerf{:,IS} = PxProdt;
%                     PxProdPerfCFPast = PxProdPerfCF{1,end};
%                     if PosDatePast > 0
%                         PxProdPerfCFPast = PxProdPerfCFPast(:,1:PosDatePast);
%                     else
%                         PxProdPerfCFPast = [];
%                     end
%                     PxProdPerfCF{1,IS} = [PxProdPerfCFPast,PxProdt];
                end
%                figure
%                plot(s(:,1),y,'.r',s(:,1),PxProdt1,'.g')
%                grid on
            end
        end
        PxProdRisk = [];
    end
    
