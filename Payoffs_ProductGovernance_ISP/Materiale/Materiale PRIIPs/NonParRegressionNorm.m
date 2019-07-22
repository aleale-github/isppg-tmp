function Fit = NonParRegressionNorm(InX,InY,OutX,delta)
% Non Parametric Regression Using Uniform Kernel
%
% Performs the Non parametric Regression Using the Kernel Method.
% The Kernel used is the Uniform's.
% The explicit formula is:
%
%  OutY(i) = ( Sum.j( K(OutX(i), InX(j), delta) )^-1 * Sum.j( InY(j) * K(OutX(i), InX(j), delta) )
%
% hence:
%
%  OutY(i) =  Sum.j( InY(j) * w(j) )
%
% where:
%
%  w(j) = K(OutX(i), InX(j), delta) / ( Sum.j( K(OutX(i), InX(j), delta) )

%     NumPoint = length(InX);
% 
% % Bucket Simulation
%     NumBucket = NumPoint/1000;
%     if NumBucket < 1
%         NumBucket = 1;
%     end
%     b = cumsum(ones(NumBucket,1),1)*1000;
%     a = [1;b(1:end-1)+1];
%     for i = 1:NumBucket
%         Ker = KernelNorm(OutX(1,a(i):b(i)),InX(1,a(i):b(i)),delta);
%         den = sum(Ker);
%         OutY = (InY(1,a(i):b(i))*Ker)./den;
%         Fit(a(i):b(i),1) = OutY;
%     end

    m = length(OutX);
    if m > 0
        sOutX = sort(OutX);
        if m > 1000
            step = (sOutX(end)-sOutX(1))/m;
            newOutX = sOutX(1):step:sOutX(end);
            Ker = KernelNorm(newOutX,InX,delta);
            den = sum(Ker);
            newOutY = (InY*Ker)./den;
            OutY = interp1(newOutX,newOutY,OutX);
        else
            Ker = KernelNorm(OutX,InX,delta);
            den = sum(Ker);
            OutY = (InY*Ker)./den;
        end
    end

    Fit = OutY;



