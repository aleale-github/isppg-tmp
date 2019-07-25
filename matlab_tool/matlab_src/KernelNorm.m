function KernelNorm = KernelNorm(X,Z,delta)
% Normal Kernel
%
% Evaluate the Uniform Kernel in the point x for all the observation Z's and with bandwidth delta 
    lX = length(X);
    lZ = length(Z);

    X = (ones(lZ,1)*X);
    Z = (ones(lX,1)*Z)';

    expon = exp(-(X-Z).^2/(2*delta^2));

    KernelNorm = 1/sqrt(pi*2*delta)*expon; 	
