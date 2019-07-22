function [RendNetFee,InvNetFee,RendCost,InvCost] = PerformanceCostScenariosPRIIPs(PriceDC,PxProdPerf,PxProdPerfCF,CostPerf,CleanPrice,ParaProd,Inv,T,FixingDate,PaymentDate,PaymentDateProd,PaymentDateInt,PosDateInt,PosDateFixing,DistrFee,ManagFee,DeltaFundingExitFee,ExitFee,Perc,Cat,Method,ProductType,Scope)

% Setting
  format long

% Product Price Scenarios
    if Cat == 2
    % Internal Variable
        NumIntScen = size(PxProdPerf,2);
        NumPerfScen = size(Perc,1);

    % Space Dimension
        RendNetFee = zeros(NumPerfScen,NumIntScen);
        InvNetFee = zeros(NumPerfScen,NumIntScen);
      
    % Check for Stock Price or Issue Price
        PxProd0 = repmat(CleanPrice+UpFrontFee,NumPerfScen,1);
        for i = 1:NumIntScen
            PxProdPerfPerc = PxProdPerf{1,i};            
            if i == NumIntScen
                Coupon = [-PxProd0,PxProdPerfPerc(:,i)*(1-T(i)*DeltaFundingExitFee)];                
            elseif i < NumIntScen
                Coupon = [-PxProd0,PxProdPerfPerc(:,i)*(1-T(i)*DeltaFundingExitFee-ExitFee)];
            end            
            for j = 1:NumPerfScen
                RendNetFee(j,i) = xirr(Coupon(j,:)',[PaymentDate(1);PaymentDate(i+1)],0.01,1000,3);
                InvNetFee(j,i) = RendNetFee(j,i).*Inv;                
            end
        end
    elseif Cat == 3
    % Internal Variable
        NumIntScen = size(PxProdPerf,2);
        NumPerfScen = size(Perc,1);
        PxProdPerfPerc = zeros(NumPerfScen,NumIntScen);

    % Space Dimension
        RendNetFee = zeros(NumPerfScen,NumIntScen);
        InvNetFee = zeros(NumPerfScen,NumIntScen);
        RendCost = zeros(3,NumIntScen);
        InvCost = zeros(1,NumIntScen);
    
    % Calculation for residual life
%         ProdLifeIHP = zeros(2,1);
%         ProdLife = yearfrac(PaymentDate(1),PaymentDate(end),12);
%         if (ProdLife-floor(ProdLife))<= 0.001
%            ProdLife = floor(ProdLife);
%         end
%         if (ProdLife > 1) && (ProdLife < 3)
%             ProdLifeIHP(1) = ProdLife-1;
%             ProdLifeIHP(2) = [];
%         elseif (ProdLife >= 3)
%             ProdLifeIHP(1) = ProdLife-1;
%             ProdLifeIHP(2) = ProdLife-ceil(ProdLife/2);
%         else
%             ProdLifeIHP(2) = [];
%             ProdLifeIHP(1) = [];
%         end  
        
        
ProdLifeIHP = zeros(2,1);
        ProdLife = yearfrac(PaymentDate(1),PaymentDate(end),12);
        if (ProdLife-floor(ProdLife))<= 0.001
           ProdLife = floor(ProdLife);
        end
        if (ProdLife > 1) && (ProdLife < 3)
            ProdLifeIHP(1) = yearfrac(PaymentDateInt(1),PaymentDate(end),2);
            ProdLifeIHP(2) = [];
        elseif (ProdLife >= 3)
            ProdLifeIHP(1) = yearfrac(PaymentDateInt(1),PaymentDate(end),2);
            ProdLifeIHP(2) = yearfrac(PaymentDateInt(2),PaymentDate(end),2);
        else
            ProdLifeIHP(2) = [];
            ProdLifeIHP(1) = [];
        end  
        
        
    % Check for Stock Price or Issue Price
        PxProd0 = repmat(CleanPrice+DistrFee,NumPerfScen,1);
        PxProd0P = repmat(CleanPrice+DistrFee,10000,1);
        %PxProd0 = repmat(CleanPrice,NumPerfScen,1);
        for i = 1:NumIntScen
        % Percentile Selection
            RefCol = size(PxProdPerf{1,i},2);
            if strcmp(Scope,'Stress')
                if i < NumIntScen
                    %if PosDateInt(i) >= PosDateFixing(end)
                    %    if T(i+1) <= 1
                    %        Perc = 1;
                    %    else
                    %        Perc = 5;
                    %    end
                    %else
                        if T(i) <= 1
                            Perc = 1;
                        else
                            Perc = 5;
                        end
                    %end
                else
                    if T(i) <= 1
                        Perc = 1;
                    else
                        Perc = 5;
                    end
                end
            end
             if strcmp(ProductType,'FRC')
              PxProdPerf{1,i}=PxProdPerf{1,i}.*PriceDC;
             end
            DFCostPerf = CostPerf{1,end};
            
            if strcmp(ProductType,'ATC1') || strcmp(ProductType,'ATC2') || strcmp(ProductType,'ATD1') || strcmp(ProductType,'ATD2') || strcmp(ProductType,'WOATD')|| strcmp(ProductType,'WOATD1')|| strcmp(ProductType,'ATCTW')
                if i == NumIntScen
                % IRR Net RHP
                    NumLifeP = size(PxProdPerfCF{:,i},2);
                    PxProdPerfCF1P = PxProdPerfCF{:,i};
                    PointerAutP = sum(PxProdPerfCF1P,2) < 1;
                % Check for Fix payment Date    
                    NumPaymentDateP = size(PaymentDate,1);
                    CouponFixP = ParaProd(end-NumPaymentDateP+1:end,1);                    
                    if CouponFixP(i) > 0
                        PointerAutP = 1;
                    end                    
                    CouponP = [-PxProd0P,PxProdPerf{1,i}-sum(ManagFee.*DFCostPerf(end-NumLifeP+1:end,1),1).*PointerAutP];  % Metodo 1
                    CouponP(:,2) = max(CouponP(:,2),0);
                elseif i < NumIntScen
                % IRR Net INT
                    ExitFeeNewP = ExitFee;
                    TotLifeP = size(PxProdPerfCF{:,end},2);                 
                    ResLifeP = TotLifeP-size(PxProdPerfCF{:,i},2);
                    NumLifeP = size(PxProdPerfCF{:,i},2);                  
                    PxProdPerfCF1P = PxProdPerfCF{:,i};
                    if isempty(PxProdPerfCF1P)
                        PointerAutP = 1;
                    else
                        PointerAutP = sum(PxProdPerfCF1P,2) < 1;
                    end
                    CouponP = [-PxProd0P,PxProdPerf{:,i}-ExitFeeNewP.*PointerAutP-DeltaFundingExitFee.*ProdLifeIHP(i).*PointerAutP-sum(ManagFee.*DFCostPerf(end-NumLifeP+1:end,1),1).*PointerAutP];  % Metodo 1
                    CouponP(:,2) = max(CouponP(:,2),0);
                end
                PxProdPerfP = CouponP(:,end)./abs(CouponP(:,1));
                else
                PxProdPerfP = PxProdPerf{1,i};    
            end
            [Move,PosScenPerc] = Percentile(PxProdPerfP,Perc,RefCol);
            Perform = PxProdPerf{1,i};
            Move = Perform(PosScenPerc,:);
            if strcmp(ProductType,'FRC')
              PxProdPerf{1,i}=PxProdPerf{1,i}./PriceDC; 
              PxProd0=PxProd0./PriceDC(PosScenPerc);
              if strcmp(Scope,'Stress')
              DistrFee=DistrFee./PriceDC(PosScenPerc);
              CleanPrice=CleanPrice./PriceDC(PosScenPerc);
              else
              DistrFee=DistrFee./PriceDC(PosScenPerc(2));
              CleanPrice=CleanPrice./PriceDC(PosScenPerc(2));
              end
            end
            
            Move(isnan(Move) == 1) = 0;
            if strcmp(ProductType,'FRC')
            Move=Move./PriceDC(PosScenPerc);
            end
            PxProdPerfPerc(:,i) = Move;
            
            if i == NumIntScen
            % IRR Net RHP
                NumLife = size(PxProdPerfCF{:,i},2);
                PxProdPerfCF1 = PxProdPerfCF{:,i};
                if strcmp(ProductType,'ATC1') || strcmp(ProductType,'ATC2') || strcmp(ProductType,'ATD1') || strcmp(ProductType,'ATD2') || strcmp(ProductType,'WOATD')|| strcmp(ProductType,'WOATD1')|| strcmp(ProductType,'ATCTW')
                %    PointerAut1 = PxProdPerfCF1 >= 1;
                %    PointerAut1 = sum(PointerAut1,2)*10000;
                %    PointerAut1(PointerAut1 == 0) = 1;
                %    PointerAut1(PointerAut1 > 100) = 0;
                %    PointerAut = PointerAut1(PosScenPerc,1);
                    PointerAut1 = sum(PxProdPerfCF1,2) < 1;
                    PointerAut = PointerAut1(PosScenPerc,1);                    
                % Check for Fix payment Date    
                    NumPaymentDate = size(PaymentDate,1);
                    CouponFix = ParaProd(end-NumPaymentDate+1:end,1);                    
                    if CouponFix(i) > 0
                        PointerAut = 1;
                    end                    
                else
                    PointerAut = 1;
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
                    strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')|| ...
                    strcmp(ProductType,'CPCA2') || strcmp(ProductType,'CPCV2') || ...
                    strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE')
                    %Coupon = [-PxProd0,PxProdPerfPerc(:,i)-ManagFee*DFCostPerf(1,end)];  % Metodo 1
                    PxProd0 = repmat(1/(1+ManagFee),NumPerfScen,1);
                    PxProd1 = ones(NumPerfScen,1);                    
                    Coupon = [-PxProd0,PxProd0.*(1+PxProdPerfPerc(:,i))];  % Metodo MOD
                    CouponNet1 = [-PxProd1,PxProd1.*(1+PxProdPerfPerc(:,i))];
                else                
                    Coupon = [-PxProd0,PxProdPerfPerc(:,i)-sum(ManagFee.*DFCostPerf(end-NumLife+1:end,1),1).*PointerAut];  % Metodo 1
                    Coupon(:,2) = max(Coupon(:,2),0);
                end
            elseif i < NumIntScen
            % IRR Net INT
                if strcmp(ProductType,'TC') || strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPF') || strcmp(ProductType,'TCS')
                    T = round(T,1);
                    NumPaymentDate = size(FixingDate,1);
                    CouponFreq = NumPaymentDate/T(end);
                     A = ParaProd(end-NumPaymentDate+1:end,1);                     
                    ExitFeeNew = A.*1/CouponFreq.*ExitFee;
                    ExitFeeNew = sum(ExitFeeNew(end-(T(end)-T(i))*CouponFreq+1:end),1);
                elseif  strcmp(ProductType,'TCA') || strcmp(ProductType,'TPPUA') || strcmp(ProductType,'TPPFA') || strcmp(ProductType,'TCSA')
                    T = round(T,1);
                    NumPaymentDate = size(FixingDate,1);
                    CouponFreq = NumPaymentDate/T(end);
                    A = ParaProd(end-NumPaymentDate+1:end,1); 
                    ExitFeeNew = A.*1/CouponFreq.*ExitFee;
                    ExitFeeNew = sum(ExitFeeNew(end-(T(end)-T(i))*CouponFreq+1:end),1);
                else
                    ExitFeeNew = ExitFee;
                end
                TotLife = size(PxProdPerfCF{:,end},2);                 
                ResLife = TotLife-size(PxProdPerfCF{:,i},2);
                NumLife = size(PxProdPerfCF{:,i},2);                  
                PxProdPerfCF1 = PxProdPerfCF{:,i};
                if strcmp(ProductType,'ATC1') || strcmp(ProductType,'ATC2') || strcmp(ProductType,'ATD1') || strcmp(ProductType,'ATD2') || strcmp(ProductType,'WOATD')||strcmp(ProductType,'WOATD1')|| strcmp(ProductType,'ATCTW')
                    %PointerAut1 = PxProdPerfCF1(:,1:end-1) >= 1;
                    %PointerAut1 = sum(PointerAut1,2)*10000;
                    %PointerAut1(PointerAut1 == 0) = 1;
                    %PointerAut1(PointerAut1 > 100) = 0;
                    %PointerAut = PointerAut1(PosScenPerc,1);
                    if isempty(PxProdPerfCF1)
                        PointerAut = 1;
                    else
                        PointerAut1 = sum(PxProdPerfCF1,2) < 1;
                        PointerAut = PointerAut1(PosScenPerc,1);
                    end
                % Check for Fix payment Date    
                    %NumPaymentDate = size(PaymentDate,1);
                    %CouponFix = ParaProd(end-NumPaymentDate+1:end,1);                    
                    %if CouponFix(i) > 0
                    %    PointerAut = 1;
                    %end
                else
                    PointerAut = 1;
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
                    strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')|| ...
                    strcmp(ProductType,'CPCA2') || strcmp(ProductType,'CPCV2') || ...
                    strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE')
                    %Coupon = [-PxProd0,PxProdPerfPerc(:,i)-ExitFeeNew-ManagFee*CostPerf{1,i}];  % Metodo 1
                    PxProd0 = repmat(1/(1+ManagFee),NumPerfScen,1);
                    PxProd1 = ones(NumPerfScen,1);
                    Coupon = [-PxProd0,PxProd0.*(1+PxProdPerfPerc(:,i)-ExitFeeNew)];  % Metodo MOD    
                    CouponNet1 = [-PxProd1,PxProd1.*(1+PxProdPerfPerc(:,i))];
                else
                    if strcmp(ProductType,'FRC')
                    ExitFeeNew=ExitFeeNew./PriceDC(PosScenPerc);
                    end
                    Coupon = [-PxProd0,PxProdPerfPerc(:,i)-ExitFeeNew.*PointerAut-DeltaFundingExitFee.*ProdLifeIHP(i).*PointerAut-sum(ManagFee.*DFCostPerf(end-NumLife+1:end,1),1).*PointerAut];  % Metodo 1
                    Coupon(:,2) = max(Coupon(:,2),0);
                end
            end
        % Value for Perfomance Table
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
                strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE')
                
                %Coupon(:,1) = Coupon(:,1)*Inv-Inv;
                %Coupon(:,2:end) = Coupon(:,2:end).*Inv;   
                %InvNetFee(:,i) = Inv.*(Coupon(:,end)./abs(Coupon(:,1))-1);
                Coupon(:,1) = -1;
                CouponNet1 = Inv.*CouponNet1;
                InvNetFee(:,i) = Inv.*(Coupon(:,end)./abs(Coupon(:,1))-1);                
            else
                InvNetFee(:,i) = Inv.*Coupon(:,end)./abs(Coupon(:,1));
            end
            for j = 1:size(Coupon(:,1),1)
                RIY = xirr(Coupon(j,:)',[PaymentDate(1);PaymentDate(i+1)],-0.9,1000,3);
                if isnan(RIY) == 1
                    RIY = (Coupon(j,end)./abs(Coupon(j,1))-1)./T(i);
                end
                RendNetFee(j,i) = RIY;   
            end

        % Value for Cost Table
            a = 2;
            if size(PosScenPerc,1) == 1
                a = 1;
            end
            if strcmp(Method,'DRM1') || strcmp(Method,'DRM3') || strcmp(Method,'DRM2')

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
                    strcmp(ProductType,'CPCKII') || strcmp(ProductType,'CPCKIE')           
                    %CouponGross = [-0,PxProdPerfPerc(a,i)];
                    %CouponGross(isnan(CouponGross) == 1) = 0;                    
                    %CouponGross(:,1) = (CouponGross(:,1)-1).*Inv;
                    %CouponGross(:,2:end) = CouponGross(:,2:end).*Inv;
                    %CouponGross(:,2:end) = (CouponGross(:,2:end)+1).*Inv;
                    %Coupon(:,1) = (PxProd0+DistrFee+ManagFee-1).*Inv;
                    %RIYGross = xirr(CouponGross',[PaymentDate(1);PaymentDate(i+1)],-0.9,1000,3);
                    CouponNet = [-Inv,Inv.*Coupon(a,2)];
                    CouponGross = CouponNet1(a,:);
                    RIYGross = xirr(CouponGross',[PaymentDate(1);PaymentDate(i+1)],-0.9,1000,3);
                    RIYNet = xirr(CouponNet',[PaymentDate(1);PaymentDate(i+1)],-0.9,1000,3);
                    RendCost(1,i) = RIYGross-RIYNet;
                    RendCost(2,i) = RIYGross;
                    RendCost(3,i) = RIYNet;
                    RendCost(isnan(RendCost(:,i)) == 1,i) = 0;                
                    if i < NumIntScen
                        InvCost(1,i) = Inv.*(DistrFee+ManagFee+ExitFeeNew);
                    else
                       InvCost(1,i) = Inv.*(DistrFee+ManagFee);
                    end                    
                else
                    CouponGross = [-CleanPrice,PxProdPerfPerc(a,i)];
                    CouponNet = Coupon(a,:);
                    RIYGross = xirr(CouponGross',[PaymentDate(1);PaymentDate(i+1)],-0.9,1000,3);
                    RIYNet = xirr(CouponNet',[PaymentDate(1);PaymentDate(i+1)],-0.9,1000,3);
                    RendCost(1,i) = RIYGross-RIYNet;
                    RendCost(2,i) = RIYGross;
                    RendCost(3,i) = RIYNet;
                    RendCost(isnan(RendCost(:,i)) == 1,i) = 0;
                    if size(PointerAut,1) > 1
                        PointerAut = PointerAut(2);
                    end
                    if i < NumIntScen
                        if strcmp(ProductType,'CWC') || strcmp(ProductType,'CWP')
                            Sp = (ParaProd(end-1)-ParaProd(end-2))/2;
                            InvCost(1,i) = Sp.*2.*Inv./ParaProd(end-1);    
                        else
                            if strcmp(ProductType,'FRC') && strcmp(Scope,'Perf')
                            ExitFeeNew=ExitFeeNew(2);
                            end
                            InvCost(1,i) = Inv./abs(CouponNet(1)).*(DistrFee+ExitFeeNew.*PointerAut+DeltaFundingExitFee.*ProdLifeIHP(i).*PointerAut+sum(ManagFee.*DFCostPerf(end-NumLife+1:end,1),1).*PointerAut);
                        end
                    else
                        if strcmp(ProductType,'CWC') || strcmp(ProductType,'CWP')
                            Sp = (ParaProd(end-1)-ParaProd(end-2))/2;
                            InvCost(1,i) = Sp.*Inv./ParaProd(end-1);
                        else
                            InvCost(1,i) = Inv./abs(CouponNet(1)).*(DistrFee+sum(ManagFee.*DFCostPerf(end-NumLife+1:end,1),1).*PointerAut);
                        end
                        
                    end
                end
            end
             if strcmp(ProductType,'FRC')
              PxProd0=PxProd0.*PriceDC(PosScenPerc);
              if strcmp(Scope,'Stress')
              DistrFee=DistrFee.*PriceDC(PosScenPerc);
              CleanPrice=CleanPrice.*PriceDC(PosScenPerc);
              else
              DistrFee=DistrFee.*PriceDC(PosScenPerc(2));
              CleanPrice=CleanPrice.*PriceDC(PosScenPerc(2));
              end    
            end
        end
    end                    
    
 