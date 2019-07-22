function [SRI,MRM,CRM,VEV] = RiskIndicatorPRIIPs(PxProd0,PxProdRisk,T,DistrFee,Rating,Cat)

% VaR Calculation
    if Cat == 2
        VaR = PxProdRisk;
    elseif Cat == 3
        VaR = prctile(PxProdRisk,2.5)./(PxProd0+DistrFee);
        %VaR = prctile(PxProdRisk,2.5);
    end
% VEV Calculation
    VEV = (sqrt(3.842-2*log(VaR))-1.96)/sqrt(T(end));
% MRM Calculation
    if VEV <= 0
        MRM = 1;
    elseif VEV < 0.005
        MRM = 1;
    elseif (VEV >= 0.005) && (VEV < 0.05)
        MRM = 2;
    elseif (VEV >= 0.05) && (VEV < 0.12)
        MRM = 3;                
    elseif (VEV >= 0.12) && (VEV < 0.20)
        MRM = 4;        
    elseif (VEV >= 0.20) && (VEV < 0.30)
        MRM = 5;        
    elseif (VEV >= 0.30) && (VEV < 0.80)
        MRM = 6;        
    elseif (VEV >= 0.80) || (isnan(VEV) == 1)
        MRM = 7;
    end
    
% CRM Calculation
% Step 1 - Credit Quality step
    if strcmp(Rating{1},'AAA') || strcmp(Rating{1},'Aaa') || strcmp(Rating{1},'AAA-')
        CRM = 0;
    elseif strcmp(Rating{1},'AA') || strcmp(Rating{1},'Aa') || strcmp(Rating{1},'Aa1') || strcmp(Rating{1},'Aa2') || strcmp(Rating{1},'Aa3') || strcmp(Rating{1},'AA-') || strcmp(Rating{1},'AA+')
        CRM = 1;
    elseif strcmp(Rating{1},'A') || strcmp(Rating{1},'A1') || strcmp(Rating{1},'A2') || strcmp(Rating{1},'A3') || strcmp(Rating{1},'A-') || strcmp(Rating{1},'A+')
        CRM = 2;
    elseif strcmp(Rating{1},'BBB') || strcmp(Rating{1},'Baa') || strcmp(Rating{1},'Baa1') || strcmp(Rating{1},'Baa2') || strcmp(Rating{1},'Baa3') || strcmp(Rating{1},'BBB-') || strcmp(Rating{1},'BBB+')
        CRM = 3;        
    elseif strcmp(Rating{1},'BB') || strcmp(Rating{1},'Ba') || strcmp(Rating{1},'Ba1') || strcmp(Rating{1},'Ba2') || strcmp(Rating{1},'Ba3') || strcmp(Rating{1},'BB-') || strcmp(Rating{1},'BB+')
        CRM = 4;        
    elseif strcmp(Rating{1},'B') || strcmp(Rating{1},'B1') || strcmp(Rating{1},'B2') || strcmp(Rating{1},'B3') || strcmp(Rating{1},'B-') || strcmp(Rating{1},'B+')
        CRM = 5;        
    else
        CRM = 6;
    end
% Step 2 - Correction for Maturity
    if T <= 1
        if (CRM >= 2) && (CRM <= 5) 
            CRM = CRM-1;
        end
    elseif (T > 12)
        if (CRM >= 4) && (CRM <= 5) 
            CRM = CRM+1;
        end
    end

% Step 3 - conversion to CRM
    if CRM == 0
        CRM = CRM+1;
    end
    
% SRI Calculation
    if CRM <= 2
        SRI = MRM;
    elseif CRM == 3
        if MRM <= 3
            SRI = 3;
        elseif MRM > 3
            SRI = MRM;
        end
    elseif (CRM == 4) && (CRM == 5)
        if MRM <= 5
            SRI = 5;
        elseif MRM > 5
            SRI = MRM;
        end
    elseif CRM == 6
        if MRM <= 6
            SRI = 6;
        elseif MRM > 6
            SRI = MRM;
        end
    end