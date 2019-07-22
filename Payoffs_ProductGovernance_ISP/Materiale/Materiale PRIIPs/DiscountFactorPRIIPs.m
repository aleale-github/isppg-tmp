
function [r0T,P0T,rtTT] = DiscountFactorPRIIPs(IRC,Today,DateFut,Fwd)  

% Check for DateFut
    if isempty(DateFut)
        r0T = [];
        P0T = [];
        rtTT = [];
    else
    % Internal Variable
        NumSpotRateDates = size(DateFut,1);

    % Interest Rate Curve (uploading)
        TenorDates = IRC(:,1);
        r0 = IRC(:,2);
        irdc = IRDataCurve('Zero',Today,TenorDates,r0','InterpMethod','linear','Compounding',-1,'Basis',12);
    % Spot Interpolated Discount Calculation    
        P0T = getDiscountFactors(irdc,DateFut);
        r0T = getZeroRates(irdc,DateFut);
        rtTT = getForwardRates(irdc,DateFut);

    % Forward Discount Factor (Drift)
        if Fwd > 0
            %NumYear = floor(Fwd/12);
            %NumMonth = Fwd-NumYear*12;
            for i = 1:NumSpotRateDates                
                FwdDatesRate = datenum(year(DateFut),month(DateFut)+Fwd,day(DateFut));
                FwdDatesRate = DateMatWrkDays(m2xdate(FwdDatesRate));
                fwd = getForwardRates(irdc,[DateFut(i);FwdDatesRate(i)]);
                rtTT(i,1) = fwd(2);
            end
        end
    end
        
    