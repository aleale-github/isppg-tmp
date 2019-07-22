function [PRIIPsNetFeeRend,PRIIPsNetFeeInv] = ScenariosPRIIPs(N,PaymentDate,PaymentDateInt,PxProd,Px0,Prc,Fee,ExitCost)

% Scenarios Generation (Money vs Return)
    
% Internal Variable
    NumSim = size(PxProd,1);
    NumSimInt = size(PxProd,2);
    Px0GrossFee = ones(NumSim,1)*(Px0+Fee+ExitCost);
    CashFlowGrossFee = [-Px0GrossFee,PxProd];
    IRRNetFee = zeros(NumSim,NumSimInt);
    
    for s = 1:NumSim
        for t = 1:NumSimInt
            if t == NumSimInt
                IRRNetFee(s,t) = xirr([CashFlowGrossFee(s,1),CashFlowGrossFee(s,1+t)]',PaymentDate,[],[],3);
            else
                PaymentDateNew = PaymentDate(PaymentDate(:,1) < PaymentDateInt(t),:);
                PaymentDateNew = [PaymentDateNew;PaymentDateInt(t)];
                IRRNetFee(s,t) = xirr([CashFlowGrossFee(s,1),CashFlowGrossFee(s,1+t)]',PaymentDateNew,[],[],3);
            end
        end
    end
    PRIIPsNetFeeRend = prctile(IRRNetFee,Prc,1);
    PRIIPsNetFeeInv = PRIIPsNetFeeRend*N;   