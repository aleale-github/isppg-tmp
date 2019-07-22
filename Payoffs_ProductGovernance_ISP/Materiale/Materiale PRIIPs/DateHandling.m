    function [FixingDateFut,IssueDate,PaymentFloatDateFut,PaymentFixDateFut,PaymentDateInt,FixingUndPast,PaymentDateIssue,PosDataNumHistSim] = DateHandling(Data,FixingDate,PaymentFloatDate,PaymentFixDate,PricingDate,SampleSize,ProductType)

% Issue Date    
    
    IssueDate = PaymentFloatDate(1);
% Correction for Pricing post Issue Date
    FixingDateFut = FixingDate(FixingDate > PricingDate);
    PaymentFloatDateFut = PaymentFloatDate(PaymentFloatDate > PricingDate);
    PaymentFixDateFut = PaymentFixDate(PaymentFixDate > PricingDate);
% Correction for Fixing before Issue Date
    FixingDatePast = FixingDate(FixingDate <= PricingDate);

% Internal variable
    NumPastDate = length(FixingDatePast);
    NumPastDate1 = PaymentFixDate(PaymentFixDate <= PricingDate);
%     NumAsset = size(Data,2)-1;
    NumAsset = size(Data,2)/2;
% Memory Space    
    FixingUndPast = zeros(NumPastDate,NumAsset);

    
% Duration in year fraction (ACT/ACT)
    PaymentDateIssue = sort([PaymentFloatDate;PaymentFixDate]);
    PaymentDateIssue(diff(PaymentDateIssue) == 0) = [];
    PaymentDateIssue(1) = [];
    PaymentDate = PaymentFixDateFut;
    if isempty(NumPastDate1) == 1
        if strcmp(ProductType,'TC') || strcmp(ProductType,'TCA') || ... 
            strcmp(ProductType,'TCS') || strcmp(ProductType,'TCSA') || ...
            strcmp(ProductType,'TPPF') || strcmp(ProductType,'TPPFA') || ...
            strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPUA') || ... 
            strcmp(ProductType,'TPPVI') || ...
            strcmp(ProductType,'CPPVI') || strcmp(ProductType,'CPPVE') || ...
            strcmp(ProductType,'CCTI') || strcmp(ProductType,'CCTE') || ...
            strcmp(ProductType,'CPFKII') || strcmp(ProductType,'CPFKIE') || ...
            strcmp(ProductType,'CPCI') || strcmp(ProductType,'CPCE') || ...
            strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')
            YMat = round(yearfrac(PaymentDate(1),PaymentDate(end),12),1);
%            if (YMat-floor(YMat))<= 5/365
%                YMat = floor(YMat);
%            end
        else
            YMat = yearfrac(PaymentDate(1),PaymentDate(end),12);
            if (YMat-floor(YMat))<= 0.001
               YMat = floor(YMat);
            end
        end
    % Date for Intermediate Scenarios (RTS Annex II. PART I, 19-20-21)
        if (YMat > 1) && (YMat < 3)
            Date1Y = datenum(year(PaymentDate(1))+1,month(PaymentDate(1)),day(PaymentDate(1)));
            Date1Y = DateMatWrkDays(m2xdate(Date1Y));
            DateIntY = [];
        elseif (YMat >= 3)
            Date1Y = datenum(year(PaymentDate(1))+1,month(PaymentDate(1)),day(PaymentDate(1)));
            DateIntY = datenum(year(PaymentDate(1))+ceil(YMat/2),month(PaymentDate(1)),day(PaymentDate(1)));
        else
            Date1Y = [];
            DateIntY = [];
        end        
    else
        %YMat = round(yearfrac(PricingDate,PaymentDate(end),12),2);
        if strcmp(ProductType,'TC') || strcmp(ProductType,'TCA') || ... 
            strcmp(ProductType,'TCS') || strcmp(ProductType,'TCSA') || ...
            strcmp(ProductType,'TPPF') || strcmp(ProductType,'TPPFA') || ...
            strcmp(ProductType,'TPPU') || strcmp(ProductType,'TPPUA') || ... 
            strcmp(ProductType,'TPPVI') || ...
            strcmp(ProductType,'CPPVI') || strcmp(ProductType,'CPPVE') || ...
            strcmp(ProductType,'CCTI') || strcmp(ProductType,'CCTE') || ...
            strcmp(ProductType,'CPFKII') || strcmp(ProductType,'CPFKIE') || ...
            strcmp(ProductType,'CPCI') || strcmp(ProductType,'CPCE') || ...
            strcmp(ProductType,'CPCA') || strcmp(ProductType,'CPCV')
            YMat = round(yearfrac(PricingDate,PaymentDate(end),12),1);
        else
            YMat = yearfrac(PricingDate,PaymentDate(end),3);
        end
    % Date for Intermediate Scenarios (RTS Annex II. PART I, 19-20-21)
        if (YMat > 1) && (YMat < 3)
            Date1Y = datenum(year(PricingDate(1))+1,month(PricingDate(1)),day(PricingDate(1)));
            DateIntY = [];
        elseif (YMat >= 3)
            Date1Y = datenum(year(PricingDate(1))+1,month(PricingDate(1)),day(PricingDate(1)));
            DateIntY = datenum(year(PricingDate(1))+ceil(YMat/2),month(PricingDate(1)),day(PricingDate(1)));
        else
            Date1Y = [];
            DateIntY = [];
        end        
    end

    PaymentDateInt = [busdate(busdate(Date1Y,1),-1);busdate(busdate(DateIntY,1),-1)]; 

% Position Date for xY simulation
    NumHistSim = size(SampleSize,1);
    PosDataNumHistSim = zeros(NumHistSim,1);
    MinDataSim = datenum(year(PricingDate(1))-SampleSize(1),month(PricingDate(1)),day(PricingDate(1)));
    for j = 1:NumAsset
        PosDataNumHistSim(j,1) = sum(Data(:,j*2-1)>= MinDataSim,1);
        
% Searching for Past Value Date
    for t = 1:NumPastDate
        if isempty(find(Data(:,j*2-1) == FixingDatePast(t))) == 0
        FixingUndPast(t,j) = Data(FixingDatePast(t) == Data(:,j*2-1),j*2);
        c = 1;
        else
        waitfor(msgbox('Non cè la data passata'));
        end
        while isempty(FixingUndPast(t,1)) == 1
            FixingUndPast = Data(FixingDatePast(t)+c == Data(:,1),2:end);
            c = c+1;
        end
    end
    end

    end
    



