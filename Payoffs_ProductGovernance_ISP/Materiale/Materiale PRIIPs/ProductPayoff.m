function [Payoff,PxUndInt,CouponFix,PriceDC] = ProductPayoff(ParaProd,PxUnd,CouponFreq,PosDateFixing,PosDateInt,NumCouponFix,FixingUndPast,FixingType,PathDep,ProductType,DataType)
% Product Payoff and Performance Scenarios Generation


% Internal Variable
    NumFixingType = size(FixingType,1);
    NumPosDateInt = size(PosDateInt,1);
    NumSim = size(PxUnd,1);
    NumAsset = size(PxUnd,2);
    FsSup = cumsum(FixingType,1);
    FsInf = [1;FsSup(1:end-1,1)+1];   
% Memory Space
    NumCycle = NumSim/100;
    B = linspace(100,NumSim,NumCycle)';
    A = [1;B(1:end-1)+1];
    PxUndInt = zeros(NumSim,NumAsset,NumPosDateInt);
    RiskFactor = zeros(100,NumAsset,NumFixingType);    
    PriceDC=ones(NumSim,1);
    for n = 1:NumCycle
        
    % Asian Mean
        c = 1;
        PxUndFixing = PxUnd(A(n):B(n),:,PosDateFixing);
        if FixingUndPast ~= 0
            FixingUndPast1 = permute(repmat(FixingUndPast',[1,1,100]),[3,1,2]);
            PxUndFixing = cat(3,FixingUndPast1,PxUndFixing);
        end

        for i = 1:NumFixingType
            RiskFactor(:,:,i) = mean(PxUndFixing(:,:,FsInf(c):FsSup(c)),3);
            c = c+1;
        end
    % Underlying for Regression
        if strcmp(PathDep{1},'D') || strcmp(PathDep{1},'W') || strcmp(PathDep{1},'M')
            if strcmp(DataType{1},'Rate')
                PxUndInt(A(n):B(n),:,:) = PxUnd(A(n):B(n),:,PosDateInt);
            else
                %PxUndInt(A(n):B(n),:,:) = PxUnd(A(n):B(n),:,PosDateInt)./repmat(RiskFactor(:,:,1),[1,1,NumPosDateInt]);
                PxUndInt(A(n):B(n),:,:) = PxUnd(A(n):B(n),:,PosDateInt);
            end
        end  
        
    
    % Product Payoff (Equity Protection Cap)
        if strcmp(ProductType{1},'MLC')
            [Payoff1,CouponFix1] = MaxLongCapPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC1')
            [Payoff1,CouponFix1] = MaxLongCap1Payoff(ParaProd,RiskFactor,NumCouponFix,100);  
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC2')
            [Payoff1,CouponFix1] = MaxLongCap2Payoff(ParaProd,RiskFactor,NumCouponFix,100);  
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC3')
            [Payoff1,CouponFix1] = MaxLongCap3Payoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC4')
            [Payoff1,CouponFix1] = MaxLongCap4Payoff(ParaProd,RiskFactor,NumCouponFix,100);   
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC5')
            [Payoff1,CouponFix1] = MaxLongCap5Payoff(ParaProd,RiskFactor,NumCouponFix,100);  
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC6')
            [Payoff1,CouponFix1] = MaxLongCap6Payoff(ParaProd,RiskFactor,NumCouponFix,100);    
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC7')
            [Payoff1,CouponFix1] = MaxLongCap7Payoff(ParaProd,RiskFactor,NumCouponFix,100); 
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'MLC8')
            [Payoff1,CouponFix1] = MaxLongCap8Payoff(ParaProd,RiskFactor,NumCouponFix,100); 
    % Product Payoff (Equity Protection Sigma Cap)
        elseif strcmp(ProductType{1},'MLCS')
            [Payoff1,CouponFix1] = MaxLongCapPayoffSigma(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff (Equity Protection Cap)
        elseif strcmp(ProductType{1},'EPS1')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = EquityProtectionShort1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);        
    % Product Payoff (Long Benchmark)
        elseif strcmp(ProductType{1},'LB')
            [Payoff1,CouponFix1] = LongBenchmarkPayoff(ParaProd,RiskFactor,NumCouponFix,100);              
    % Product Payoff (Bonus Cap)
        elseif strcmp(ProductType{1},'BC1')
            BarrierFactor = PxUnd(A(n):B(n),:,end);        
            [Payoff1,CouponFix1] = BonusCap1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff (Bonus Cap)
        elseif strcmp(ProductType{1},'BC2')
            BarrierFactor = PxUnd(A(n):B(n),:,end);        
            [Payoff1,CouponFix1] = BonusCap2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff (Bonus Cap)
        elseif strcmp(ProductType{1},'BC3')
            BarrierFactor = PxUnd(A(n):B(n),:,end);        
            [Payoff1,CouponFix1] = BonusCap3Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);             
    % Product Payoff (Fix to Floater)        
        elseif strcmp(ProductType{1},'F2F')
            [Payoff1,CouponFix1] = FixToFloaterPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);
    % Product Payoff (Fix to Floater 1)        
        elseif strcmp(ProductType{1},'F2F1')
            [Payoff1,CouponFix1] = FixToFloaterPayoff1(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);
 % Product Payoff (Fix to Floater 2)        
        elseif strcmp(ProductType{1},'F2F2')
            [Payoff1,CouponFix1] = ReverseFloaterCapped(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);
    % Product Payoff (Fix to Floater 1)        
        elseif strcmp(ProductType{1},'IL')
            [Payoff1,CouponFix1] = InflationLinkedPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100); 
    % Product Payoff (Inflation Linked Note)        
        elseif strcmp(ProductType{1},'IL1')
            [Payoff1,CouponFix1] = InflationLinkedPayoff1(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);         
    % Product Payoff (Gap Long Certificate)      
        elseif strcmp(ProductType{1},'GAP')
            [Payoff1,CouponFix1] = GapLongCertificatePayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);     
    % Product Payoff (Worst Off)        
        elseif strcmp(ProductType{1},'WO')
            [Payoff1,CouponFix1] = WorstOffPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff (Worst Off)        
        elseif strcmp(ProductType{1},'WO1')
            [Payoff1,CouponFix1] = WorstOff1Payoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff (Worst Off 2, Digital Memory)        
        elseif strcmp(ProductType{1},'WO2')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = WorstOff2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff (Worst Off Autocallable)        
        elseif strcmp(ProductType{1},'WOATD')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = WorstOffAutocallablePayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff (Worst Off Autocallable)        
        elseif strcmp(ProductType{1},'WOATD1')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = WorstOffAutocallablePayoff1(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff (TYPE 1 AUTOCALLABLE)        
        elseif strcmp(ProductType{1},'ATC1')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = Autocallable1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff (TYPE 2 AUTOCALLABLE)        
        elseif strcmp(ProductType{1},'ATC2')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = Autocallable2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,100);            
    % Product Payoff (TYPE 1 AUTOCALLABLE DIGITAL)        
        elseif strcmp(ProductType{1},'ATD1')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AutocallableDigital1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,100);       
    % Product Payoff (TYPE 2 AUTOCALLABLE DIGITAL)        
        elseif strcmp(ProductType{1},'ATD2')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AutocallableDigital2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,100);   
     % Product Payoff (TYPE 3 AUTOCALLABLE DIGITAL)        
        elseif strcmp(ProductType{1},'ATD3')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AutocallableDigital3Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,100);    
    % Product Payoff (TYPE 1 Autocallable Twin Win)        
        elseif strcmp(ProductType{1},'ATCTW')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AutocallableTwinWinPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,100);    
    % Product Payoff (IRS Fix Payer)
        elseif strcmp(ProductType{1},'TC')
            [Payoff1,CouponFix1] = TassoCertoAmmoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);            
    % Product Payoff (IRS Fix Payer Amortizing)
        elseif strcmp(ProductType{1},'TCA')
            [Payoff1,CouponFix1] = TassoCertoAmmoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);
    % Product Payoff (IRS Tasso Certo Super)
        elseif strcmp(ProductType{1},'TCS')
            [Payoff1,CouponFix1] = TassoCertoSuperPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);            
    % Product Payoff (IRS Tasso Certo Super Ammo)
        elseif strcmp(ProductType{1},'TCSA')
            [Payoff1,CouponFix1] = TassoCertoSuperPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);
    % Product Payoff (Tasso Protetto Premio Unico)    
        elseif strcmp(ProductType{1},'TPPU')
            [Payoff1,CouponFix1] = TassoProtettoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);
    % Product Payoff (Tasso Protetto Premio Unico Ammo)
        elseif strcmp(ProductType{1},'TPPUA')
            [Payoff1,CouponFix1] = TassoProtettoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);            
    % Product Payoff (Tasso Protetto Premio Frazionato)
        elseif strcmp(ProductType{1},'TPPF')
            [Payoff1,CouponFix1] = TassoProtettoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);  
    % Product Payoff (Tasso Protetto Premio Frazionato Ammo)
        elseif strcmp(ProductType{1},'TPPFA')
            [Payoff1,CouponFix1] = TassoProtettoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);            
    % Product Payoff (Cambio Protetto Plain Vanilla Import)        
       elseif strcmp(ProductType{1},'CPPVI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);   
    % Product Payoff (Cambio Protetto Plain Vanilla Export)        
       elseif strcmp(ProductType{1},'CPPVE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Protetto Collar Import)        
       elseif strcmp(ProductType{1},'CPCI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);   
    % Product Payoff (Cambio Protetto Collar Export)        
       elseif strcmp(ProductType{1},'CPCE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
     % Product Payoff (Cambio Protetto Collar Knock In Import)        
       elseif strcmp(ProductType{1},'CPCKII')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);   
    % Product Payoff (Cambio Protetto Collar Knock In Export)        
       elseif strcmp(ProductType{1},'CPCKIE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Protetto Forward Knock In Import)        
       elseif strcmp(ProductType{1},'CPFKII')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Protetto Forward Knock In Export)        
       elseif strcmp(ProductType{1},'CPFKIE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Protetto Best Forward Import)        
       elseif strcmp(ProductType{1},'CPBFI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);   
    % Product Payoff (Cambio Protetto Best Forward Export)        
       elseif strcmp(ProductType{1},'CPBFE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);     
    % Product Payoff (Cambio Protetto Forward Knock In Knock Out Import)        
       elseif strcmp(ProductType{1},'CPFKIKOI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Protetto Forward Knock In Knock Out Export)        
       elseif strcmp(ProductType{1},'CPFKIKOE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);              
    % Product Payoff (Cambio Certo Termine Import)        
       elseif strcmp(ProductType{1},'CCTI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);  
    % Product Payoff (Cambio Certo Termine Export)        
       elseif strcmp(ProductType{1},'CCTE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Certo Termine Import)        
       elseif strcmp(ProductType{1},'CCPI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);  
    % Product Payoff (Cambio Certo Termine Export)        
       elseif strcmp(ProductType{1},'CCPE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);  
    % Product Payoff (Prezzo Certo Termine Import - EURUSD)        
       elseif strcmp(ProductType{1},'CPCA')
            Type = 'A';
            [Payoff1,CouponFix1] = PrezzoCertoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Prezzo Certo Termine Export - EURUD)        
       elseif strcmp(ProductType{1},'CPCV')
            Type = 'V';
            [Payoff1,CouponFix1] = PrezzoCertoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Prezzo Certo Termine Import - USDEUR)        
       elseif strcmp(ProductType{1},'CPCA2')
            Type = 'A';
            [Payoff1,CouponFix1] = PrezzoCertoPayoff2(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Prezzo Certo Termine Export - USDEUR)        
       elseif strcmp(ProductType{1},'CPCV2')
            Type = 'V';
            [Payoff1,CouponFix1] = PrezzoCertoPayoff2(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Dinamico Accumulator Import)        
       elseif strcmp(ProductType{1},'CDAI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioDinamicoAccumulatorPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Dinamico Accumulator Import)        
       elseif strcmp(ProductType{1},'CDAE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioDinamicoAccumulatorPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Dinamico Accumulator Leva Import)        
       elseif strcmp(ProductType{1},'CDALI')
            Type = 'I';
            [Payoff1,CouponFix1] = CambioDinamicoAccumulatorPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);
    % Product Payoff (Cambio Dinamico Accumulator Leva Export)        
       elseif strcmp(ProductType{1},'CDALE')
            Type = 'E';
            [Payoff1,CouponFix1] = CambioDinamicoAccumulatorPayoff(ParaProd,RiskFactor,NumCouponFix,Type,100);              
    % Product Payoff Annual Digital Tipo 1        
        elseif strcmp(ProductType{1},'AD1')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AnnualDigital1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff Annual Digital  Tipo 2
        elseif strcmp(ProductType{1},'AD2')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AnnualDigital2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff Annual Digital  Tipo 3
        elseif strcmp(ProductType{1},'AD3')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AnnualDigital3Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);            
    % Product Payoff Annual Digital  Tipo 4
        elseif strcmp(ProductType{1},'AD4')
            [Payoff1,CouponFix1] = AnnualDigital4Payoff(ParaProd,RiskFactor,NumCouponFix,100);   
    % Product Payoff Annual Digital  Tipo 5
        elseif strcmp(ProductType{1},'AD5')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AnnualDigital5Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff Annual Digital  Tipo 6
        elseif strcmp(ProductType{1},'AD6')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = AnnualDigital6Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff Annual Digital  Tipo 7
        elseif strcmp(ProductType{1},'AD7')
            BarrierFactor = PxUnd(A(n):B(n),:,end);              
            [Payoff1,CouponFix1] = AnnualDigital7Payoff(ParaProd,RiskFactor,NumCouponFix,100);             
    % Product Payoff Annual Digital  Tipo 8
        elseif strcmp(ProductType{1},'AD8')
            BarrierFactor = PxUnd(A(n):B(n),:,end);              
            [Payoff1,CouponFix1] = AnnualDigital8Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);   
    % Product Payoff Annual Digital  Tipo 9
        elseif strcmp(ProductType{1},'AD9')
            BarrierFactor = PxUnd(A(n):B(n),:,end);              
            [Payoff1,CouponFix1] = AnnualDigital9Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
    % Product Payoff Annual Digital Short            
        elseif strcmp(ProductType{1},'ADS')
            BarrierFactor = PxUnd(A(n):B(n),:,end);            
            [Payoff1,CouponFix1] = AnnualDigitalShortPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);                
    % Product Payoff Double Barrier Digital
        elseif strcmp(ProductType{1},'DBD')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = BarrierDigitalPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);  
    % Product Payoff Double Barrier Digital Tipo 2
        elseif strcmp(ProductType{1},'DBD2')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = BarrierDigital2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);  
    % Product Payoff Max Long Cap Digital
        elseif strcmp(ProductType{1},'MLCD')
            BarrierFactor = PxUnd(A(n):B(n),:,end);
            [Payoff1,CouponFix1] = MaxLongCapDigitalPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100); 
    % Product Payoff Twin Win
        elseif strcmp(ProductType{1},'TW')
            [Payoff1,CouponFix1] = TwinWinPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff Twin Win Knock-Out
        elseif strcmp(ProductType{1},'TWKO')
            [Payoff1,CouponFix1] = TwinWinKnockOutPayoff(ParaProd,RiskFactor,NumCouponFix,100);        
    % Product Payoff Accelerator
        elseif strcmp(ProductType{1},'A')
            [Payoff1,CouponFix1] = AcceleratorPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff Accelerator Tipo 3
        elseif strcmp(ProductType{1},'A3')
            [Payoff1,CouponFix1] = Accelerator3Payoff(ParaProd,RiskFactor,NumCouponFix,100);        
    % Product Payoff Covered Warrant Call
        elseif strcmp(ProductType{1},'CWC')
            [Payoff1,CouponFix1] = CoveredWarrantCallPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff Covered Warrant Put
        elseif strcmp(ProductType{1},'CWP')
            [Payoff1,CouponFix1] = CoveredWarrantPutPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff Covered Warrant Put
        elseif strcmp(ProductType{1},'OR')
            [Payoff1,CouponFix1] = ObiettivoRendimentoPayoff(ParaProd,RiskFactor,NumCouponFix,100);             
     % Product Payoff (Swap Note)        
        elseif strcmp(ProductType{1},'CMS')
            [Payoff1,CouponFix1] = ConstantMaturitySwapNotePayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);  
     % Product Payoff Podium Equity            
        elseif strcmp(ProductType{1},'PEL')
            [Payoff1,CouponFix1] = PodiumEquityLinkedPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,100);     
    % Product Fix Rate Currency
        elseif strcmp(ProductType{1},'FRC')
            [Payoff1,CouponFix1,PriceDC1] = FixRateEURXXXPayoff(ParaProd,RiskFactor,NumCouponFix,100);
    % Product Payoff (LECOIP)
        elseif strcmp(ProductType{1},'LECOIP')
            BarrierFactor = PxUnd(A(n):B(n),:,end);        
            [Payoff1,CouponFix1] = LECOIPPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,100);
     % Product Payoff Obbligazione Strutturata-Basket Currencies Linked         
        elseif strcmp(ProductType{1},'OS1')
            [Payoff1,CouponFix1] = ObbligazioneStrutturata1Payoff(ParaProd,RiskFactor,NumCouponFix,100);            
        end            
        Payoff(A(n):B(n),:) = Payoff1;
        CouponFix(A(n):B(n),:) = CouponFix1;
        if strcmp(ProductType{1},'FRC')
        PriceDC(A(n):B(n),:)=PriceDC1;
        end
    end
    
    
% ------------------------------------------------------------------------------------------------------------------------------------------------------------------- %    
   
    function [Payoff,CouponFix] = MaxLongCapPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = squeeze(min(max(ParaProd(3).*(RF-1)+1,ParaProd(2)),ParaProd(3).*(ParaProd(1)-1)+1));

    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = MaxLongCap1Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
        
    % ParaProd(1) : Percentuale Cap
    % ParaProd(2) : Percentuale di Protezione
    % ParaProd(3) : Fattore di Partecipazione (su performance positiva)
    % ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix) : Weights
    % ParaProd(end-NumCouponFix+1:end,1) : Cedole Fisse

    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = squeeze(min(max(ParaProd(3).*(RF-1)+1,ParaProd(2)),ParaProd(1)));

    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
    
    function [Payoff,CouponFix] = MaxLongCap2Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = squeeze(min(ParaProd(1),ParaProd(2)+ParaProd(3).*max(0,RF-ParaProd(2))));
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);    
        
    function [Payoff,CouponFix] = MaxLongCap3Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = squeeze(min(ParaProd(1),max(ParaProd(2),RF+(ParaProd(3)-1).*max(RF-1,0))));
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = MaxLongCap4Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = squeeze(min(max(ParaProd(3).*(1-RF)+1,ParaProd(2)),ParaProd(3).*(ParaProd(1)-1)+1));

    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);        
        
    function [Payoff,CouponFix] = MaxLongCap5Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = squeeze(sum(RF.*Weights,2));
        CheckKI = sum(RF >= ParaProd(3),2);
        CheckKI(CheckKI > 0) = 1;
        LP = ParaProd(2)+CheckKI.*(ParaProd(3)-ParaProd(2));
        
        Payoff = squeeze(min(max(ParaProd(4).*(RF(:,end)-1)+1,LP),ParaProd(4).*(ParaProd(1)-1)+1));
       

    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);  
        
    function [Payoff,CouponFix] = MaxLongCap6Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = ParaProd(2) + max(RF-ParaProd(2),0) + (ParaProd(3)-1).*max(RF-ParaProd(1),0);
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
  
    function [Payoff,CouponFix] = MaxLongCap7Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = RF;
        Payoff(RF <= ParaProd(2),end) = ParaProd(2);
        Payoff(RF > ParaProd(2) & RF < ParaProd(4)) = RF(RF > ParaProd(2) & RF < ParaProd(4))+(ParaProd(3)-1).*(RF(RF > ParaProd(2) & RF < ParaProd(4))-1);
        Payoff(RF >= ParaProd(4),end) = ParaProd(1);
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);    
        
    function [Payoff,CouponFix] = MaxLongCap8Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
    % Index Calculation
        VRI = min(RiskFactor(:,:,1:end-NumCouponFix),[],3);
        RF = squeeze(RiskFactor(:,:,end)./VRI);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',size(RF));
        RF = sum(RF.*Weights,2);
        Payoff = squeeze(min(ParaProd(1),max(ParaProd(2),RF+(ParaProd(3)-1).*max(RF-1,0))));
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);

    function [Payoff,CouponFix] = MaxLongCapPayoffSigma(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = RF+ParaProd(5);
        Payoff(RF >= ParaProd(2),end) = min(ParaProd(1),max(ParaProd(3),RF(RF>ParaProd(2))+(ParaProd(4)-1).*max((RF(RF>ParaProd(2))-1),0)));
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
    
        function [Payoff,CouponFix] = EquityProtectionShort1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim)
        
        % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
        NumCouponFloat = size(RF,2);
    % Check Digital
        CheckDigital = PayoffRF <= repmat(ParaProd(1:NumCouponFloat)',NumSim,1);
        Payoff = CheckDigital.*repmat(ParaProd(NumCouponFloat+1:NumCouponFloat*2)',NumSim,1);
    % Maturity    
        IP=ParaProd(NumCouponFloat*2+1);
        Prot=ParaProd(NumCouponFloat*2+2);

        Payoff(BF > ParaProd(NumCouponFloat*2+3),end) = Payoff(BF > ParaProd(NumCouponFloat*2+3),end) + max(Prot,2-BF(BF > ParaProd(NumCouponFloat*2+3),end));
        Payoff(BF <= ParaProd(NumCouponFloat*2+3),end) = Payoff(BF <= ParaProd(NumCouponFloat*2+3),end) + IP;

     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); 
        
    function [Payoff,CouponFix] = LongBenchmarkPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff = RF;
          
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = BonusCap1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % ParaProd(1) : Livello Barriera
    % ParaProd(2) : Coupon a scadenza
    % ParaProd(end-NumCouponFix+1:end,1) : Cedole Fisse
        
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]);
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));        
        Payoff = RF;
        Payoff(RF >= ParaProd(1),end) = ParaProd(2)+1;
         
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = BonusCap2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]);
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));        
        Payoff = RF;
        Payoff(BF >= ParaProd(2),end) = max(Payoff(BF >= ParaProd(2),end),ParaProd(3)+1);
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = BonusCap3Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]);
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));        
        Payoff = RF;
        Payoff(BF >= ParaProd(1),end) = min(ParaProd(3),max(Payoff(BF >= ParaProd(1),end),ParaProd(2)));
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);           
        
    function [Payoff,CouponFix] = FixToFloaterPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    % Internal Variable
        RF = squeeze(RiskFactor);
        NumCoupon = size(RF,2);
    % Deal Parameters adjstment    
        Spread = repmat(ParaProd(1:NumCoupon,1)',NumSim,1); % Spread
        PF = repmat(ParaProd(NumCoupon+1:NumCoupon*2,1)',NumSim,1); % Partecipation Factor
        Floor = repmat(ParaProd(NumCoupon*2+1:NumCoupon*3,1)',NumSim,1); % Floor
        Cap = repmat(ParaProd(NumCoupon*3+1:NumCoupon*4,1)',NumSim,1); % Cap
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFreq = repmat(CouponFreq',NumSim,1);

    % Fix Coupon    
        CouponFix = CpnFix.*CouponFreq;
        
    % Payoff
        Payoff = min(max(PF.*RF+Spread,Floor),Cap);
        
        if NumCoupon < NumCouponFix
            Pos = ParaProd(end-NumCouponFix+1:end,1) > 0;
            CouponFreq(:,Pos) = [];
        end       
        Payoff = Payoff.*CouponFreq;
        Payoff(:,end) = Payoff(:,end)+1;
        
    function [Payoff,CouponFix] = FixToFloaterPayoff1(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
 
    % Internal Variable
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);    
        RF = squeeze(RF);
        NumCoupon = size(RF,2);
    % Deal Parameters adjstment    
        Spread = repmat(ParaProd(1:NumCoupon,1)',NumSim,1); % Spread
        PF = repmat(ParaProd(NumCoupon+1:NumCoupon*2,1)',NumSim,1); % Partecipation Factor
        Floor = repmat(ParaProd(NumCoupon*2+1:NumCoupon*3,1)',NumSim,1); % Floor
        Cap = repmat(ParaProd(NumCoupon*3+1:NumCoupon*4,1)',NumSim,1); % Cap
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFreq = round(repmat(CouponFreq',NumSim,1));
    
    % Fix Coupon    
        CouponFix = CpnFix.*CouponFreq;
        
    % Payoff
        Payoff = min(max(PF.*(RF-1)+Spread,Floor),Cap);
        Payoff = Payoff.*CouponFreq(:,end-NumCoupon+1:end);
        Payoff(:,end) = Payoff(:,end)+1;
        
function [Payoff,CouponFix] = ReverseFloaterCapped (ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    % Internal Variable
        RF = squeeze(RiskFactor);
        NumCoupon = size(RF,2);
    % Deal Parameters adjstment    
        Diff = repmat(ParaProd(1:NumCoupon,1)',NumSim,1); % Spread
        PF = repmat(ParaProd(NumCoupon+1:NumCoupon*2,1)',NumSim,1); % Partecipation Factor
        Floor = repmat(ParaProd(NumCoupon*2+1:NumCoupon*3,1)',NumSim,1); % Floor
        Cap = repmat(ParaProd(NumCoupon*3+1:NumCoupon*4,1)',NumSim,1); % Cap
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFreq = repmat(CouponFreq',NumSim,1);

    % Fix Coupon    
        CouponFix = CpnFix.*CouponFreq;

    % Payoff

        Payoff = min(Diff-PF.*RF*(365/360),Cap);
        
        if NumCoupon < NumCouponFix
            Pos = ParaProd(end-NumCouponFix+1:end,1) > 0;
            CouponFreq(:,Pos) = [];
        end       
        Payoff = Payoff.*CouponFreq;
        Payoff(:,end) = Payoff(:,end)+1;
        
    function [Payoff,CouponFix] = InflationLinkedPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    % Internal Variable
        RF = squeeze(RiskFactor(:,:,2:end)./RiskFactor(:,:,1:end-1));
        NumCoupon = size(RF,2);
        
    % Deal Parameters adjstment    
        PF = repmat(ParaProd(1),NumSim,NumCoupon); % Partecipation Factor
        Floor = repmat(ParaProd(2),NumSim,NumCoupon); % Floor
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFreq = repmat(CouponFreq',NumSim,1);

    % Fix Coupon    
        CouponFix = CpnFix; 
    % Payoff
        Payoff = PF.*max(Floor,RF-1);       
        Payoff(:,end) = Payoff(:,end)+1;
        
     function [Payoff,CouponFix] = InflationLinkedPayoff1(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    % Internal Variable
        RF = squeeze(RiskFactor(:,:,2:end)./RiskFactor(:,:,1:end-1));
        NumCoupon = size(RF,2);
        
    % Deal Parameters adjstment    
        PF = repmat(ParaProd(1),NumSim,NumCoupon); % Partecipation Factor
        Floor = repmat(ParaProd(2),NumSim,NumCoupon); % Floor
        Cap = repmat(ParaProd(3),NumSim,NumCoupon); %Cap
        Spread=repmat(ParaProd(4),NumSim,NumCoupon); %Spread
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFreq = repmat(CouponFreq',NumSim,1);

    % Fix Coupon    
        CouponFix = CpnFix; 
    % Payoff
        Payoff = min(Cap,PF.*max(Floor,(RF-1)+Spread));       
        Payoff(:,end) = Payoff(:,end)+1;

    function [Payoff,CouponFix] = GapLongCertificatePayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    
    % Index Calculation    
        RF = RiskFactor(:,:,2:end)./RiskFactor(:,:,1:end-1);
        RF = squeeze(RF);
        
    % Gap Calculation
        Barrier = repmat(ParaProd(1),size(RF));
        BarrierEvent = RF < Barrier;
        GapIndicator = cumsum(BarrierEvent,2);
        GapIndicator = BarrierEvent.*GapIndicator;
        GapIndicator(GapIndicator>1) = 0;
        Gap = RF.*GapIndicator;
        GapEnd = sum(Gap,2);
        GapIndicatorEnd = sum(GapIndicator,2);
                    
        
    % Payoff
        Payoff = (1-GapIndicatorEnd)*ParaProd(3) + GapIndicatorEnd.*max(0, 1 + ParaProd(2)*(GapEnd-ParaProd(1)));
        
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);  

    function [Payoff,CouponFix] = WorstOffPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)    
          
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]);
        Payoff1 = sum(RF >= ParaProd(2),2);
        WO = min(RF,[],2);
        Payoff(Payoff1 == 2,:) = ParaProd(1);
        Payoff(Payoff1 < 2,:) = WO(Payoff1 < 2,:);
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = WorstOff1Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)    
    

    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]);
        NumDigEvent = size(RF,3);
        Payoff = zeros(NumSim,NumDigEvent);
        for t = 1:NumDigEvent
            WO = min(RF(:,:,t),[],2);
            if t < NumDigEvent
                Payoff(WO >= ParaProd(t),t) = ParaProd(NumDigEvent+1+t);
                Payoff(WO < ParaProd(t),t) = 0;
            elseif t == NumDigEvent
                Payoff(WO >= ParaProd(t),t) = ParaProd(NumDigEvent+1+t)+ParaProd(NumDigEvent+2+t);
                Payoff(WO < ParaProd(t),t) = max(WO(WO < ParaProd(t)),ParaProd(NumDigEvent+3+t));            
            end
        end        
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);        

    function [Payoff,CouponFix] = WorstOff2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
       
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        %PayoffRF = RF;
        WO = squeeze(min(RF,[],2));
        WOB = min(BF,[],2);
        
    % Digital Effect
        CouponDigital = repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        DigitalEffect = (WO >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        
    % Memory Effect
        PointerMemoryEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix)',NumSim,1);
        CouponMemory = repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        MF = (WO >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponMemory;
        MemoryEffect = zeros(100,size(PointerMemoryEffect,2));
        MemoryEffect1 = zeros(100,size(PointerMemoryEffect,2));
        MemoryEffect1(MF(:,1) == 0,1) = CouponMemory(MF(:,1) == 0,1);
        for j = 2:size(PointerMemoryEffect,2)
            MemoryEffect1(MF(:,j) == 0,j) = MemoryEffect1(MF(:,j) == 0,j-1)+CouponMemory(MF(:,j) == 0,j);
            %MemoryEffect(MF(:,j-1) < MF(:,j),j) = MemoryEffect1(MF(:,j-1) < MF(:,j),j-1)+CouponMemory(MF(:,j-1) < MF(:,j),j);
            MemoryEffect(MF(:,j-1) < MF(:,j),j) = MemoryEffect1(MF(:,j-1) < MF(:,j),j-1);
        end
        MemoryEffect = MemoryEffect.*PointerMemoryEffect;
        %MemoryEffect = MemoryEffect;
        
     % Final Coupon T
        %PointerAutoCallEnd = sum(PointerAutoCall(:,1:end-1),2);
        %PointerAutoCallEnd = PointerAutoCall(:,end);
        %PayoffEnd = ones(size(BF,1),1).*ParaProd(2+NumAutoCall)+1;
        Payoff = DigitalEffect+MemoryEffect;
        PayoffEnd = DigitalEffect(:,end)+MemoryEffect(:,end)+1;
        
        Payoff(:,end) = PayoffEnd;
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = WorstOffAutocallablePayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim)    

    %{
        ParaProd(1:NumCouponFix) = Barriere Digital
        ParaProd(NumCouponFix+1:NumCouponFix*2) = Barriere Autocall
        ParaProd(NumCouponFix*2+1:NumCouponFix*3) = Coupon Digital
        ParaProd(NumCouponFix*3+1) = Rimborso a scadenza sopra la Barriera
            Autocall
        ParaProd(NumCouponFix*3+2) = Rimborso a scadenza sopra la Barriera
            Digital
        ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2) = Date effetto
            Digital
        ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1) = Date di
            Autocall
        ParaProd(end-NumCouponFix+1:end,1) = Cedole Fix
    %}
        
    % Index Calculation
        RF1 = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        
        NumEvent = size(RF1,3);        
    % Correction for Worst of
        RF = zeros(NumSim,NumCouponFix);
        for t = 1:NumEvent
            WO = min(RF1(:,:,t),[],2);
            RF(:,t) = WO;
            if t == NumEvent
                BF1 = squeeze(BarrierFactor./RiskFactor(:,:,1));                
                WO1 = min(BF1,[],2);
                BF = WO1;
            end
        end

    % Digital Effect
        PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponDigital = repmat(ParaProd(NumCouponFix*2+1:NumCouponFix*3)',NumSim,1);
        DigitalEffect = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        DigitalEffect = DigitalEffect.*PointerDigitalEffect;
        %DigitalEffect = DigitalEffect.*CouponFreq;  
        
    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = DigitalEffect;
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall == 0) = 0;
    % AutoCall Coupon T
%        PointerAutoCallEnd = cumsum(PointerAutoCall,2);
        PointerAutoCallEnd = PointerAutoCall(:,end);
        PayoffEnd = BF;
%        PayoffEnd(BF < ParaProd(NumCouponFix*2+1)) = BF(BF < ParaProd(NumCouponFix*2+1),end);
        PayoffEnd(BF >= ParaProd(NumCouponFix*2)) = ParaProd(NumCouponFix*3+1);
        PayoffEnd(BF >= ParaProd(NumCouponFix)) = ParaProd(NumCouponFix*3+2);        
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0;
        AutoCallCoupon(:,end) = PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;        
   
    function [Payoff,CouponFix] = WorstOffAutocallablePayoff1(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim)    
 
    % Index Calculation
        RF1 = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        
        NumEvent = size(RF1,3);        
    % Correction for Worst of
        RF = zeros(NumSim,NumCouponFix);
        for t = 1:NumEvent
            WO = min(RF1(:,:,t),[],2);
            RF(:,t) = WO;
            if t == NumEvent
                BF1 = squeeze(BarrierFactor./RiskFactor(:,:,1));                
                WO1 = min(BF1,[],2);
                BF = WO1;
            end
        end

    % Digital Effect
        PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponDigital = repmat(ParaProd(NumCouponFix*2+1:NumCouponFix*3)',NumSim,1);
        DigitalEffect = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        DigitalEffect = DigitalEffect.*PointerDigitalEffect;
        %DigitalEffect = DigitalEffect.*CouponFreq;  
     
        % Memory Effect
        PointerMemoryEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponMemory = repmat(ParaProd(NumCouponFix*2+1:NumCouponFix*3)',NumSim,1);
        MF = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponMemory;
        MemoryEffect = zeros(100,size(PointerMemoryEffect,2));
        MemoryEffect1 = zeros(100,size(PointerMemoryEffect,2));
        MemoryEffect1(MF(:,1) == 0,1) = CouponMemory(MF(:,1) == 0,1);
        for j = 2:size(PointerMemoryEffect,2)
            MemoryEffect1(MF(:,j) == 0,j) = MemoryEffect1(MF(:,j) == 0,j-1)+CouponMemory(MF(:,j) == 0,j);
            %MemoryEffect(MF(:,j-1) < MF(:,j),j) = MemoryEffect1(MF(:,j-1) < MF(:,j),j-1)+CouponMemory(MF(:,j-1) < MF(:,j),j);
            MemoryEffect(MF(:,j-1) < MF(:,j),j) = MemoryEffect1(MF(:,j-1) < MF(:,j),j-1);
        end
        MemoryEffect = MemoryEffect.*PointerMemoryEffect;    
        
        
    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = DigitalEffect+MemoryEffect;
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall == 0) = 0;
    % AutoCall Coupon T
%        PointerAutoCallEnd = cumsum(PointerAutoCall,2);
        PointerAutoCallEnd = PointerAutoCall(:,end);
        PayoffEnd = BF;
%        PayoffEnd(BF < ParaProd(NumCouponFix*2+1)) = BF(BF < ParaProd(NumCouponFix*2+1),end);
        PayoffEnd(BF >= ParaProd(NumCouponFix*2)) = ParaProd(NumCouponFix*3+1);
        PayoffEnd = PayoffEnd + AutoCallCoupon(:,end);
        %PayoffEnd(BF >= ParaProd(NumCouponFix)) = ParaProd(NumCouponFix*3+2);        
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0;
        AutoCallCoupon(:,end) = PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;  
        
    function [Payoff,CouponFix] = Autocallable1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim)

    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]);
        BF = BarrierFactor./RiskFactor(:,:,1);
    % Pointer Matrix
        PointerAutoCall = squeeze(RF >= ParaProd(1));
        
    % Internal Variable
        NumAutoCall = size(PointerAutoCall,2);

        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;        
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 0;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = repmat(ParaProd(NumAutoCall+2:end-NumCouponFix-1,1)',NumSim,1);
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall ~= 1) = 0;
    % AutoCall Coupon T
        PointerAutoCallEnd = cumsum(PointerAutoCall,2);
%         PayoffEnd = BF;
        PayoffEnd = squeeze(RF(:,:,end));
%         PayoffEnd(BF >= ParaProd(NumAutoCall+1)) = ParaProd(end-NumCouponFix)+1;
%         PayoffEnd(BF >= ParaProd(NumAutoCall)) = ParaProd(end-NumCouponFix-1)+1;
          PayoffEnd(squeeze(RF(:,:,end)) >= ParaProd(NumAutoCall+1)) = ParaProd(end-NumCouponFix)+1;
          PayoffEnd(squeeze(RF(:,:,end)) >= ParaProd(NumAutoCall)) = ParaProd(end-NumCouponFix-1)+1;

        PayoffEnd(PointerAutoCallEnd(:,end) ~= 0) = 0;
        AutoCallCoupon(:,end) = AutoCallCoupon(:,end)+PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;
  
    function [Payoff,CouponFix] = Autocallable2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,NumSim)
    %{  
        ParaProd(1:NumCouponFix) : Barriere Digital/Memoria
        ParaProd(NumCouponFix+1:NumCouponFix*2) : Digital Amounts
        ParaProd(NumCouponFix*2+1:NumCouponFix*3) : Coupon Memory
        ParaProd(NumCouponFix*3+1:NumCouponFix*4) : Barriere Autocall
        ParaProd(NumCouponFix*4+1) : Barriera a scadenza (KO)
        ParaProd(end-NumCouponFix*4+1:end-NumCouponFix*3) : Pointer
                                                            Digital Effect
        ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2) : Pointer Memory
                                                            Effect
        ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1) : Pointer
                                                            Autocall Effect
        ParaProd(end-NumCouponFix+1:end,1) : Cedole Fisse
        
    %}
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
    % CouponFrequency
        CouponFreq = repmat(CouponFreq',NumSim,1);
    % Digital Effect
        PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*4+1:end-NumCouponFix*3)',NumSim,1);
        CouponDigital = repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        DigitalEffect = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        DigitalEffect = DigitalEffect.*PointerDigitalEffect;
        %DigitalEffect = DigitalEffect.*CouponFreq;
    % Memory Effect
        PointerMemoryEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponMemory = repmat(ParaProd(NumCouponFix*2+1:NumCouponFix*3)',NumSim,1);
        MF = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponMemory;
        MemoryEffect = zeros(100,size(PointerMemoryEffect,2));
        MemoryEffect1 = zeros(100,size(PointerMemoryEffect,2));
        MemoryEffect1(MF(:,1) == 0,1) = CouponMemory(MF(:,1) == 0,1);
        for j = 2:size(PointerMemoryEffect,2)
            MemoryEffect1(MF(:,j) == 0,j) = MemoryEffect1(MF(:,j) == 0,j-1)+CouponMemory(MF(:,j) == 0,j);
            %MemoryEffect(MF(:,j-1) < MF(:,j),j) = MemoryEffect1(MF(:,j-1) < MF(:,j),j-1)+CouponMemory(MF(:,j-1) < MF(:,j),j);
            MemoryEffect(MF(:,j-1) < MF(:,j),j) = MemoryEffect1(MF(:,j-1) < MF(:,j),j-1);
        end
        MemoryEffect = MemoryEffect.*PointerMemoryEffect;
        %MemoryEffect = MemoryEffect;
        
    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(NumCouponFix*3+1:NumCouponFix*4)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;        
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;        
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
        PointerAutoCallEnd = PointerAutoCall(:,end);
    % AutoCall Coupon T-
        PointerAutoCall = PointerAutoCall.*PointerAutocallEffect;
        AutoCallCoupon = DigitalEffect+MemoryEffect;
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall == 0) = 0;
    % Final Coupon T
        %PointerAutoCallEnd = sum(PointerAutoCall(:,1:end-1),2);
        %PointerAutoCallEnd = PointerAutoCall(:,end);
        %PayoffEnd = ones(size(BF,1),1).*ParaProd(2+NumAutoCall)+1;
        PayoffEnd = DigitalEffect(:,end)+MemoryEffect(:,end)+1;
        PayoffEnd(BF < ParaProd(NumCouponFix*4+1)) = BF(BF < ParaProd(NumCouponFix*4+1),end);
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0;
        AutoCallCoupon(:,end) = AutoCallCoupon(:,end)+PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;        
        
    function [Payoff,CouponFix] = AutocallableDigital1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,NumSim)

    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));

    % CouponFrequency
        CouponFreq = repmat(CouponFreq',NumSim,1);
    % Digital Effect
        PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponDigital = repmat(ParaProd(NumCouponFix*2+2:NumCouponFix*3+1)',NumSim,1);
        DigitalEffect = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        DigitalEffect = DigitalEffect.*PointerDigitalEffect;
        %DigitalEffect = DigitalEffect.*CouponFreq;  
        
    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = DigitalEffect;
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall == 0) = 0;
    % AutoCall Coupon T
%        PointerAutoCallEnd = cumsum(PointerAutoCall,2);
        PointerAutoCallEnd = PointerAutoCall(:,end);
        PayoffEnd = DigitalEffect(:,end)+1;
        PayoffEnd(BF < ParaProd(NumCouponFix*2+1)) = BF(BF < ParaProd(NumCouponFix*2+1),end);
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0;
        AutoCallCoupon(:,end) = PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;   
        
    function [Payoff,CouponFix] = AutocallableDigital2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,NumSim)

    % Index Calculation
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
        RF1 = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF1 = squeeze(BarrierFactor./RiskFactor(:,:,1));
        Weights = repmat(ParaProd(end-NumCouponFix*3 - NumRiskFactor + 1:end-NumCouponFix*3)',[NumSim,1,NumFixing-1]);
        
        RF = squeeze(sum(RF1.*Weights,2));
        BF = squeeze(sum(BF1.*Weights(:,:,end),2));

    % CouponFrequency
        CouponFreq = repmat(CouponFreq',NumSim,1);
    % Digital Effect
        PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponDigital = repmat(ParaProd(NumCouponFix*2+2:NumCouponFix*3+1)',NumSim,1);
        DigitalEffect = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        DigitalEffect = DigitalEffect.*PointerDigitalEffect;
        %DigitalEffect = DigitalEffect.*CouponFreq;  
        
    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = DigitalEffect;
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall == 0) = 0;
    % AutoCall Coupon T
%        PointerAutoCallEnd = cumsum(PointerAutoCall,2);
        PointerAutoCallEnd = PointerAutoCall(:,end);
        PayoffEnd = DigitalEffect(:,end)+1;
        PayoffEnd(BF < ParaProd(NumCouponFix*2+1)) = BF(BF < ParaProd(NumCouponFix*2+1),end);
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0;
        AutoCallCoupon(:,end) = PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;

    function [Payoff,CouponFix] = AutocallableDigital3Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,NumSim)

    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));

    % CouponFrequency
        CouponFreq = repmat(CouponFreq',NumSim,1);
    % Digital Effect
        PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CouponDigital = repmat(ParaProd(NumCouponFix*2+3:NumCouponFix*3+2)',NumSim,1);
        CouponAutocall = repmat(ParaProd(NumCouponFix*3+3:NumCouponFix*4+2)',NumSim,1);
        DigitalEffect = (RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*CouponDigital;
        DigitalEffect = DigitalEffect.*PointerDigitalEffect;
        DigitalEffect = DigitalEffect.*(RF < repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1));
        %DigitalEffect = DigitalEffect.*CouponFreq;  
        
    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix*1)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = DigitalEffect;
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+CouponAutocall(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall == 0) = 0;
    % AutoCall Coupon T
%        PointerAutoCallEnd = cumsum(PointerAutoCall,2);
        PointerAutoCallEnd = PointerAutoCall(:,end);
        PayoffEnd = BF;                 
        PayoffEnd(BF >= ParaProd(NumAutoCall*2+1)) = ParaProd(end-NumCouponFix*3-1)+1;
        PayoffEnd(BF >= ParaProd(NumAutoCall*2+2)) = ParaProd(end-NumCouponFix*3)+1;
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0; 
        AutoCallCoupon(:,end) = AutoCallCoupon(:,end)+PayoffEnd;

    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;


    function [Payoff,CouponFix] = AutocallableTwinWinPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,CouponFreq,NumSim)

    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));

    % CouponFrequency
        CouponFreq = repmat(CouponFreq',NumSim,1);

    % Autocall Effect
        PointerAutocallEffect = repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix)',NumSim,1);        
    % Pointer Matrix
        PointerAutoCall = RF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        PointerAutoCall(PointerAutocallEffect == 0) = 0;
        NumAutoCall = size(PointerAutoCall,2);
        for a = 1:NumAutoCall
            PointerAutoCall = cumsum(PointerAutoCall,2);
        end
        PointerAutoCall(PointerAutoCall > 1) = 30000;
        PointerAutoCall(PointerAutoCall(:,end) == 0,:) = 10000;
        PointerAutoCall(PointerAutoCall == 0) = 20000;
        PointerAutoCall(PointerAutoCall == 10000) = 3;
        PointerAutoCall(PointerAutoCall == 20000) = 2;
        PointerAutoCall(PointerAutoCall == 30000) = 0;
    % AutoCall Coupon T-
        AutoCallCoupon = repmat(ParaProd(NumCouponFix+5:end-NumCouponFix*2,1)',NumSim,1); 
        AutoCallCoupon(PointerAutoCall == 1) = AutoCallCoupon(PointerAutoCall == 1)+1;
        AutoCallCoupon(PointerAutoCall ~= 1) = 0;
    % AutoCall Coupon T
        PointerAutoCallEnd = PointerAutoCall(:,end);
        PayoffEnd = RF(:,end);
        PayoffEnd(BF >= ParaProd(NumCouponFix+1)) = ParaProd(NumCouponFix+2) + max(ParaProd(NumCouponFix+3) - BF(BF >= ParaProd(NumCouponFix+1)), 0)  + max(BF(BF >= ParaProd(NumCouponFix+1)) - ParaProd(NumCouponFix+3), 0) - max(BF(BF >= ParaProd(NumCouponFix+1)) - ParaProd(NumCouponFix+4), 0);  
        % PayoffEnd(BF < ParaProd(NumCouponFix*2+1)) = BF(BF < ParaProd(NumCouponFix*2+1),end);
        PayoffEnd(PointerAutoCallEnd(:,end) ~= 3) = 0;
        AutoCallCoupon(:,end) = PayoffEnd;
    
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        Payoff = AutoCallCoupon;  
        if prod(sum(Payoff ~= 0 , 2)) ~= 1
            error('Mess with Autocall')
        end
        
    function [Payoff,CouponFix] = TassoCertoAmmoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
        
        RF = squeeze(RiskFactor);
        Payoff = (RF-repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*repmat(ParaProd(NumCouponFix+1:end,1)',NumSim,1).*repmat(CouponFreq',NumSim,1);
        Payoff(:,end) = Payoff(:,end);
        
        CouponFix = zeros(NumSim,NumCouponFix);

    function [Payoff,CouponFix] = TassoCertoSuperPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
        
        RF = squeeze(RiskFactor);
        Payoff = (max(RF,repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1))-repmat(ParaProd(1:NumCouponFix)',NumSim,1)).*repmat(ParaProd(NumCouponFix*2+1:end,1)',NumSim,1).*repmat(CouponFreq',NumSim,1);
        Payoff(:,end) = Payoff(:,end);
        
        CouponFix = zeros(NumSim,NumCouponFix);
        
    function [Payoff,CouponFix] = TassoProtettoPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
        RF = squeeze(RiskFactor);
        P = repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        K1 = repmat(ParaProd(NumCouponFix+1:NumCouponFix*2,1)',NumSim,1);
        A = repmat(ParaProd(NumCouponFix*2+1:end,1)',NumSim,1);
        dt = repmat(CouponFreq',NumSim,1);
        Payoff = (-P+max(RF-K1,0)).*A.*dt;        
       
        CouponFix = zeros(NumSim,NumCouponFix);

    function [Payoff,CouponFix] = CambioProtettoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,NumSim)
        RF = squeeze(RiskFactor);
        K1 = repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        K2 = repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        KID = ParaProd(end-5);
        KOD = ParaProd(end-4);
        KIU = ParaProd(end-3);
        KOU = ParaProd(end-2);
        L1 = ParaProd(end-1);
        L2 = ParaProd(end);
        
        if Type == 'I'
            checkKID = sum(RF(:,2:end) <= KID,2)*10000;
            checkKID(checkKID > 0) = 1;
            checkKOD = sum(RF(:,2:end) <= KOD,2)*10000;
            checkKOD(checkKOD == 0) = 1;
            checkKOD(checkKOD > 1) = 0;
            checkKIU = sum(RF(:,2:end) >= KIU,2)*10000;
            checkKIU(checkKIU > 0) = 1;            
            checkKOU = sum(RF(:,2:end) >= KOU,2)*10000;
            checkKOU(checkKOU == 0) = 1;
            checkKOU(checkKOU > 1) = 0;
            %checkKOU = 1-checkKOU;
            
            Payoff = (1/L1).*checkKID.*checkKOD.*max(K1-RF(:,end),0)-(1/L2).*checkKIU.*checkKOU.*max(RF(:,end)-K2,0);        
        elseif Type == 'E'
            checkKID = sum(RF(:,2:end) <= KID,2)*10000;
            checkKID(checkKID > 0) = 1;
            checkKOD = sum(RF(:,2:end) <= KOD,2)*10000;
            checkKOD(checkKOD == 0) = 1;
            checkKOD(checkKOD > 1) = 0;
            checkKIU = sum(RF(:,2:end) >= KIU,2)*10000;
            checkKIU(checkKIU > 0) = 1;            
            checkKOU = sum(RF(:,2:end) >= KOU,2)*10000;
            checkKOU(checkKOU == 0) = 1;
            checkKOU(checkKOU > 1) = 0;
            %checkKOU = 1-checkKOU;
            
            Payoff = (1/L2).*checkKIU.*checkKOU.*max(RF(:,end)-K2,0)-(1/L1).*checkKID.*checkKOD.*max(K1-RF(:,end),0);  
        end
        
    % Payoff in EURO       
        Payoff(:,end) = Payoff(:,end)./RF(:,end);
       
        CouponFix = zeros(NumSim,NumCouponFix);   
        
    function [Payoff,CouponFix] = CambioDinamicoAccumulatorPayoff(ParaProd,RiskFactor,NumCouponFix,Type,NumSim)
        RF = squeeze(RiskFactor);
        K1 = repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        K2 = repmat(ParaProd(NumCouponFix+1:NumCouponFix*2)',NumSim,1);
        L1 = ParaProd(end-2);
        L2 = ParaProd(end-1);
        A = ParaProd(end);
        N = size(RF,2)-1;
        if Type == 'I'
            if L1 ~= 1 || L2 ~= 1
                checkA = sum(RF(:,2:end) >= A,2);
                p = sum(RF(:,2:end) >= K1,2);
                n = checkA-p;                
                WeightP = n/N; 
                                
                c = sum(RF(:,2:end)>= K2,2);
                p = c;
                WeightC = p/N;
            else
                n = sum(RF(:,2:end) >= A,2);
                WeightP = n/N;
                WeightC = n/N;
            end
            Payoff = (1/L1).*WeightP.*max(K1-RF(:,end),0)-(1/L2).*WeightC.*max(RF(:,end)-K2,0);
            
        elseif Type == 'E'

        end
        
    % Payoff in XXX       
        Payoff(:,end) = Payoff(:,end)./RF(:,end);
       
        CouponFix = zeros(NumSim,NumCouponFix);   
        
    function [Payoff,CouponFix] = PrezzoCertoPayoff(ParaProd,RiskFactor,NumCouponFix,Type,NumSim)
        RFt = squeeze(RiskFactor(:,1,:));
        if size(RiskFactor,2) == 2
            EURXXX = squeeze(RiskFactor(:,2,:));
        else
            EURXXX = ones(size(RiskFactor,1),size(RiskFactor,3));
        end
        C1 = repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        if Type == 'A'
            Payoff = (RFt./EURXXX-C1)./(C1*NumCouponFix);
        elseif Type == 'V'
            Payoff = (C1-RFt./EURXXX)./(C1*NumCouponFix);
        end
       
        CouponFix = zeros(NumSim,NumCouponFix);
        
 function [Payoff,CouponFix] = PrezzoCertoPayoff2(ParaProd,RiskFactor,NumCouponFix,Type,NumSim)
        RFt = squeeze(RiskFactor(:,1,:));
        if size(RiskFactor,2) == 2
            EURXXX = squeeze(RiskFactor(:,2,:));
        else
            EURXXX = ones(size(RiskFactor,1),size(RiskFactor,3));
        end
        C1 = repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        if Type == 'A'
            Payoff = (RFt.*EURXXX-C1)./(C1*NumCouponFix);
        elseif Type == 'V'
            Payoff = (C1-RFt.*EURXXX)./(C1*NumCouponFix);
        end
       
        CouponFix = zeros(NumSim,NumCouponFix);
        
    function [Payoff,CouponFix,PriceDC1] = FixRateEURXXXPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
        EURXXX = squeeze(RiskFactor(:,1,2:end));
        C1 = repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        PriceDC1=RiskFactor(:,1,1);
        Payoff = C1;
        Payoff(:,end) = Payoff(:,end)+1;
        Payoff = Payoff./EURXXX;
        CouponFix = zeros(NumSim,NumCouponFix);           
        
    function [Payoff,CouponFix] = AnnualDigital1Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
        NumCouponFloat = size(RF,2);
    % Check Digital
        CheckDigital = PayoffRF >= repmat(ParaProd(1:NumCouponFloat)',NumSim,1);
        Payoff = CheckDigital.*repmat(ParaProd(end-(NumCouponFix+NumCouponFloat)+1:end-NumCouponFix)',NumSim,1);
        Payoff(BF >= ParaProd(NumCouponFloat+1),end) = 1+ParaProd(end-NumCouponFix);
        Payoff(BF < ParaProd(NumCouponFloat+1),end) = BF(BF < ParaProd(NumCouponFloat+1),end);

     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);      
        
    function [Payoff,CouponFix] = AnnualDigital2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
    
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Memory effetct
        CheckMemory = cumprod(PayoffRF >= ParaProd(1),2);
        Payoff = CheckMemory.*repmat(ParaProd(3:end-NumCouponFix)',NumSim,1);
        Payoff(BF >= ParaProd(2),end) = Payoff(BF >= ParaProd(2),end)+1;
        Payoff(BF < ParaProd(2),end) = BF(BF < ParaProd(2),end);

     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);     
        
    function [Payoff,CouponFix] = AnnualDigital3Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
       
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Memory effect
        NumCoupon = size(RF,1);
        CheckMemory = PayoffRF >= ParaProd(1);
        
        for c = 1:NumCoupon
            CheckMemory = cumsum(CheckMemory,2);
        end
        [~,j] = find(CheckMemory == 1);
        CheckMemory(CheckMemory == 1) = CheckMemory(CheckMemory == 1).*j;
        CheckMemory(CheckMemory > (NumCoupon-1)) = 1;
        Payoff = CheckMemory.*repmat(ParaProd(3:end-NumCouponFix)',NumSim,1);
        Payoff(BF >= ParaProd(1),end) = Payoff(BF >= ParaProd(1),end)+1;
        Payoff(BF < ParaProd(1),end) = Payoff(BF < ParaProd(1),end)+max(ParaProd(2),BF(BF < ParaProd(1),end));
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
    function [Payoff,CouponFix] = AnnualDigital4Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
         
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        PayoffRF = RF;
    % Digital Effect
        CheckDigital = PayoffRF >= ParaProd(1);
        Payoff = CheckDigital.*repmat(ParaProd(2:end-NumCouponFix)',NumSim,1);
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        CouponFix(:,end) = CouponFix(:,end)+1;     
        
    function [Payoff,CouponFix] = AnnualDigital5Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim)
         
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Digital Effect
        NumCoupon = size(RF,2);
        n = 1:1:NumCoupon;
        n = repmat(n,NumSim,1);
        Payoff = min(repmat(ParaProd(4:end-NumCouponFix)',NumSim,1),max(ParaProd(3),(PayoffRF-1)./n));
        
        Payoff(BF >= ParaProd(1),end) = Payoff(BF >= ParaProd(1),end)+1;
        Payoff(BF < ParaProd(1),end) = Payoff(BF < ParaProd(1),end)+max(ParaProd(2),BF(BF < ParaProd(1),end));
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        CouponFix(:,end) = CouponFix(:,end);        
        
    function [Payoff,CouponFix] = AnnualDigital7Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        PayoffRF = RF;
    % Check Digital
        CheckDigital = PayoffRF <= ParaProd(1);
        Payoff = CheckDigital.*repmat(ParaProd(2:NumCouponFix+1)',NumSim,1);
        Payoff(:,end) = Payoff(:,end)+ParaProd(end-NumCouponFix);
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);              

    function [Payoff,CouponFix] = AnnualDigital6Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Check Digital
        CheckDigital = PayoffRF >= ParaProd(1);
        Payoff = CheckDigital.*repmat(ParaProd(4:end-NumCouponFix)',NumSim,1);
        Payoff(BF >= ParaProd(2),end) = Payoff(BF >= ParaProd(2),end)+1;
        Payoff(BF < ParaProd(2),end) = max(BF(BF < ParaProd(2),end),ParaProd(3));

     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);  
        
    function [Payoff,CouponFix] = AnnualDigital8Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % ParaProd(1:NumCouponFix) : Barriere Digital
    % ParaProd(NumCouponFix+1) : Barriera a Scadenza per cedola
    % ParaProd(NumCouponFix+2) : Barriera a scadenza per sottostante
    % ParaProd(end-NumCouponFix*2+1:end-NumCouponFix) : Cedole Digital
    % ParaProd(end-NumCouponFix+1:end,1) : Cedole FIsse
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Check Digital
        CheckDigital = PayoffRF >= repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        Payoff = CheckDigital.*repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix)',NumSim,1);
%         Payoff(BF >= ParaProd(NumCouponFix+1),end) = Payoff(BF >= ParaProd(NumCouponFix+1),end)+1;
%         Payoff(BF < ParaProd(NumCouponFix+1),end) = BF(BF < ParaProd(NumCouponFix+1),end);
%         Payoff(BF < ParaProd(NumCouponFix+2),end) = ParaProd(NumCouponFix+2);

        Payoff(RF(:,end) >= ParaProd(NumCouponFix+1),end) = Payoff(BF >= ParaProd(NumCouponFix+1),end)+1;
        Payoff(RF(:,end) < ParaProd(NumCouponFix+1),end) = BF(BF < ParaProd(NumCouponFix+1),end);
        Payoff(RF(:,end) < ParaProd(NumCouponFix+2),end) = ParaProd(NumCouponFix+2);
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); 
        
    function [Payoff,CouponFix] = AnnualDigitalShortPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Check Digital

        CheckDigital = PayoffRF <= repmat(ParaProd(1:NumCouponFix)',NumSim,1);
        Payoff = CheckDigital.*repmat(ParaProd(end-NumCouponFix*2+1:end-NumCouponFix)',NumSim,1);
        Payoff(BF <= ParaProd(NumCouponFix+1),end) = Payoff(BF <= ParaProd(NumCouponFix+1),end)+1;
        Payoff(BF > ParaProd(NumCouponFix+1),end) = max(ParaProd(end-NumCouponFix*2),1-(BF(BF > ParaProd(NumCouponFix+1),1)-1));
        

     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); 

    function [Payoff,CouponFix] = BarrierDigitalPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
 
    % Internal Variable
        NumCouponFloat = size(RiskFactor,3)-1; 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Check Digital
        CheckDigital = PayoffRF >= repmat(ParaProd(1:NumCouponFloat)',NumSim,1);
        Payoff = CheckDigital.*repmat(ParaProd(NumCouponFloat+1:end-NumCouponFix)',NumSim,1);
        Payoff(BF >= ParaProd(NumCouponFloat),end) = Payoff(BF >= ParaProd(NumCouponFloat),end)+1;
        %Payoff(BF < ParaProd(NumCouponFloat+1),end) = max(BF(BF < ParaProd(NumCouponFloat+1),end),ParaProd(NumCouponFloat+2));
        Payoff(BF < ParaProd(NumCouponFloat),end) = BF(BF < ParaProd(NumCouponFloat),end);
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);                      
    
    function [Payoff,CouponFix] = BarrierDigital2Payoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
 
    % Internal Variable
        NumCouponFloat = size(RiskFactor,3)-1; 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Check Digital
        CheckDigital = PayoffRF >= repmat(ParaProd(1:NumCouponFloat)',NumSim,1);
        Payoff = CheckDigital.*repmat(ParaProd(NumCouponFloat+2:end-NumCouponFix-1)',NumSim,1);
        Payoff(BF >= ParaProd(NumCouponFloat+1),end) = Payoff(BF >= ParaProd(NumCouponFloat+1),end) + ParaProd(end-NumCouponFix);
        %Payoff(BF < ParaProd(NumCouponFloat+1),end) = max(BF(BF < ParaProd(NumCouponFloat+1),end),ParaProd(NumCouponFloat+2));
        Payoff(BF < ParaProd(NumCouponFloat+1),end) = Payoff(BF < ParaProd(NumCouponFloat+1),end) + BF(BF < ParaProd(NumCouponFloat+1),end);
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);                      
   
    function [Payoff,CouponFix] = MaxLongCapDigitalPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
        
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    % Check Digital
       % PointerDigitalEffect = repmat(ParaProd(end-NumCouponFix*3+1:end-NumCouponFix*2)',NumSim,1);
        CheckDigital = PayoffRF(:,1:end-1) <= ParaProd(1);
        CheckDigitalTot = sum(CheckDigital,2);
        CheckDigitalTot(CheckDigitalTot >= 1) = 1;
        Payoff = zeros(NumSim,NumCouponFix);
        Payoff(:,1) = CheckDigitalTot * ParaProd(end-NumCouponFix);
        Payoff(:,end) = min(ParaProd(4),max(RF(end),ParaProd(3)));
        Payoff(BF < ParaProd(2),end) = BF(BF < ParaProd(2),end);
       
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);     
        
    function [Payoff,CouponFix] = TwinWinPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
         
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        Payoff = RF(:,end);
    % Memory effetct
        NumBarrierEvent = size(RF,2);
        Payoff1 = zeros(size(RF,1),1);
        BarrierEvent = sum(RF >= ParaProd(1),2);
        CheckBarrierEvent = zeros(NumSim,1);
        CheckBarrierEvent(BarrierEvent == NumBarrierEvent,1) = 1;
        Payoff1(CheckBarrierEvent == 1,1) = Payoff(CheckBarrierEvent == 1,1);
        Payoff1(Payoff1 >= 1,1) = min((Payoff1(Payoff1 >= 1,1)-1).*ParaProd(4)+1,ParaProd(2));
        Payoff1(Payoff1 < 1,1) = -(Payoff1(Payoff1 < 1,1)-1).*ParaProd(5)+1;
        Payoff(CheckBarrierEvent == 1,1) = Payoff1(CheckBarrierEvent == 1,1);
        Payoff(CheckBarrierEvent == 0,1) = max(Payoff(CheckBarrierEvent == 0,1),ParaProd(3));

        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);

    function [Payoff,CouponFix] = TwinWinKnockOutPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    %{
        ParaProd(1): Strike Short
        ParaProd(2): Strike Long
        ParaProd(3): Floor Short
        ParaProd(4): Floor Long
        ParaProd(5): Konck Out Short
        ParaProd(6): Konck Out Long
        ParaProd(7): Partecipation Factor Short
        ParaProd(8): Partecipation Factor Long
        ParaProd(9): Rimborso a scadenza
        ParaProd(end-NumCouponFix+1:end,1): Cedole Fisse
    %}
        
    % Index Calculation
        BF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        RF = squeeze(RiskFactor(:,:,end)./RiskFactor(:,:,1));
       
        BarrierEventShort = prod(BF >= ParaProd(5),2);
        BarrierEventLong = prod(BF <= ParaProd(6),2);
        
        Payoff = ParaProd(9) + BarrierEventShort.*max(ParaProd(3), ParaProd(7).*(ParaProd(1) - RF)) + BarrierEventLong.*max(ParaProd(4), ParaProd(8).*(RF - ParaProd(2)));
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);    
        
    function [Payoff,CouponFix] = AcceleratorPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
         
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        Payoff = RF;
    % Barrier Event
        BarrierEvent = RF >= ParaProd(1);
        Payoff(BarrierEvent == 1,1) = 1+(Payoff(BarrierEvent == 1,1)-1).*ParaProd(2);
        Payoff(BarrierEvent == 0,1) = 1+(Payoff(BarrierEvent == 0,1)-1).*ParaProd(3);
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);  
        
function [Payoff,CouponFix] = Accelerator3Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
         
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        Payoff = RF;
    % Barrier Event
        BarrierEvent = RF >= ParaProd(1);
        Payoff(BarrierEvent == 1,1) = min(ParaProd(4),1+(Payoff(BarrierEvent == 1,1)-1).*ParaProd(2));
        Payoff(BarrierEvent == 0,1) = max(ParaProd(5),1+(Payoff(BarrierEvent == 0,1)-1).*ParaProd(3));
        
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
        
    function [Payoff,CouponFix] = CoveredWarrantCallPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
             
    % Index Calculation
        FRV = squeeze(RiskFactor(:,:,2:end));
        %IRV = squeeze(RiskFactor(:,:,1));
        IRV = ParaProd(3);
    % Payoff
        Payoff = max(FRV-IRV,ParaProd(1))*ParaProd(2);

    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);

    function [Payoff,CouponFix] = CoveredWarrantPutPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
             
    % Index Calculation
        FRV = squeeze(RiskFactor(:,:,2:end));
        IRV = ParaProd(3);
    % Payoff
        Payoff = max(IRV-FRV,ParaProd(1))*ParaProd(2);

    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);            

    function [Payoff,CouponFix] = ConstantMaturitySwapNotePayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    % Internal Variable
        RF = squeeze(RiskFactor);
        RF = RF/100;
        NumCoupon = size(RF,2);
    % Deal Parameters adjstment    
        PF = repmat(ParaProd(1:NumCoupon,1)',NumSim,1); % Partecipation Factor
        Floor = repmat(ParaProd(NumCoupon+1:NumCoupon*2,1)',NumSim,1); % Floor
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFreq = repmat(CouponFreq',NumSim,1);

    % Fix Coupon    
        CouponFix = CpnFix.*CouponFreq;
        
    % Payoff
        Payoff = max(PF.*RF,Floor);
        if NumCoupon < NumCouponFix
            Pos = ParaProd(end-NumCouponFix+1:end,1) > 0;
            CouponFreq(:,Pos) = [];
        end
        Payoff = Payoff.*CouponFreq;
        Payoff(:,end) = Payoff(:,end)+1;        
     
    function [Payoff,CouponFix] = ObiettivoRendimentoPayoff(ParaProd,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3)-1;
        
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing]);
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing]);
        RF = squeeze(sum(RF.*Weights,2));
        Cap = repmat(ParaProd(3:NumFixing+2)',[NumSim,1]);
        Floor = repmat(ParaProd(NumFixing+3:2*NumFixing+2)',[NumSim,1]);
        t = repmat(linspace(1,NumFixing,NumFixing),[NumSim,1]);
        AR = (RF).^(1./t)-1;
        Payoff = squeeze(min(max(AR,Floor),Cap));
        Payoff(RF(:,end) >= ParaProd(1),end) = Payoff(RF(:,end) >= ParaProd(1),end)+1;
        Payoff(RF(:,end) < ParaProd(1),end) = max(ParaProd(2),RF(RF(:,end) < ParaProd(1),end));
    % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);        
    
    function [Payoff,CouponFix] = PodiumEquityLinkedPayoff(ParaProd,CouponFreq,RiskFactor,NumCouponFix,NumSim)
    
    % Internal Variable
        NumRiskFactor = size(RiskFactor,2);
        NumFixing = size(RiskFactor,3);
    % Index Calculation
        RF = RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,NumFixing-1]);
        Bar=repmat(ParaProd(1:NumFixing-1)',[NumSim,1,NumRiskFactor]);
        Bar = permute(Bar,[1,3,2]);
        CheckBarrier = RF >= Bar;
        CheckPodium=prod(CheckBarrier,2);
        CheckPodium=squeeze(CheckPodium);
        SumNum=cumsum(ones(NumSim,NumFixing-1),2);
        Z = SumNum.*CheckPodium;
        for k = 1:size(Z,1)
            z1 = Z(k,Z(k,:)>0);
            z2 = diff(z1);
            z1(2:end) = z2;
            Z(k,Z(k,:)>0) = z1;
        end
        
    % Fix Coupon    
        CpnFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1); % Coupon Fix
        CouponFix = CpnFix.*CouponFreq';  
        
      %  if NumCoupon < NumCouponFix
      %      Pos = ParaProd(end-NumCouponFix+1:end,1) > 0;
      %      CouponFreq(:,Pos) = [];
      %  end 
        
        Payoff = Z.*ParaProd(end-NumCouponFix-NumFixing+2:end-NumCouponFix)';
        Payoff(:,end) = Payoff(:,end) +1;        
        
    function [Payoff,CouponFix] = LECOIPPayoff(ParaProd,RiskFactor,BarrierFactor,NumCouponFix,NumSim) 
  
    % Index Calculation
        RF = squeeze(RiskFactor(:,:,2:end)./repmat(RiskFactor(:,:,1),[1,1,size(RiskFactor,3)-1]));
        BF = squeeze(BarrierFactor./RiskFactor(:,:,1));
        PayoffRF = RF;
    
        CheckStrike = max(PayoffRF -1, ParaProd(1));
        Payoff = ParaProd(2) + ParaProd(3)*ParaProd(4)*mean(CheckStrike,2);
 
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);
        
function [Payoff,CouponFix] =  ObbligazioneStrutturata1Payoff(ParaProd,RiskFactor,NumCouponFix,NumSim)  %MODIFICATO   
    % Index Calculation
        NumFixing = size(RiskFactor,3); 
        NumRiskFactor = size(RiskFactor,2);
        RF=squeeze((RiskFactor(:,:,1)-RiskFactor(:,:,2:end))./RiskFactor(:,:,1));
        Weights = repmat(ParaProd(end-NumCouponFix-NumRiskFactor+1:end-NumCouponFix)',[NumSim,1,NumFixing-1]);
        RF = sum(RF.*Weights,2);
        Payoff =ParaProd(4)+ParaProd(4)*min(ParaProd(2),max(ParaProd(1),ParaProd(3)*max(0,RF)+ParaProd(5)));
  
     % Fix Coupon    
        CouponFix = repmat(ParaProd(end-NumCouponFix+1:end,1)',NumSim,1);         
        
        