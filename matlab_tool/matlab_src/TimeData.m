function [SimDateFut,PosDateFixing,PosDateInt,PosDatePayFloat,PaymentDate,PaymentDateFut,WrkDays,NumWrkDays,CouponFreq,T] = TimeData(PricingDate,PaymentFloatDateFut,PaymentFixDateFut,PaymentFloatDate,PaymentFixDate,PaymentDateInt,FixingDateFut,TimeFreq,BasisCouponFreq,ProductType)

    [PaymentDateFut,PosDatePayFloat] = DataSorting(PaymentFloatDateFut,PaymentFixDateFut);
    PaymentDate = DataSorting(PaymentFloatDate,PaymentFixDate);
    PaymentDate(PaymentDate(:,1) == 0) = [];
    CouponFreq = yearfrac(PaymentDate(1:end-1,1),PaymentDate(2:end,1),BasisCouponFreq);
    [SimDateFut,PosDateFixing,PosDateInt] = DataSorting(FixingDateFut,PaymentDateInt);
    [WrkDays,NumWrkDays] = TimeSteps(PricingDate,SimDateFut,TimeFreq,ProductType);
    %PaymentDate = [PricingDate;PaymentDateInt;PaymentDateFut(end)];
    if PricingDate < PaymentDate(1)
        PaymentDate = [PaymentDateFut(1);PaymentDateInt;PaymentDateFut(end)];
    elseif PricingDate == PaymentDate(1)
        PaymentDate = [PaymentDate(1);PaymentDateInt;PaymentDateFut(end)];
    else
        PaymentDate = [PricingDate;PaymentDateInt;PaymentDateFut(end)];
    end
   [~,~,T] = TimeSteps(PaymentDate(1),PaymentDate(2:end,1),TimeFreq,ProductType);    
%    if FixingUndPast > 0
%    else
%        PaymentDate = [PaymentDateFut(1);PaymentDateInt;PaymentDateFut(end)];
%    end
       

