nscen = 10000;


from =0;
to = 6.28*2;
step = (to-from)/(nscen-1);
X=from:step:to;

delta = .25;
Y_true  =  cos(X)+X;
Y= Y_true +randn(1,nscen);

plot(X,Y,'.y',X,NonParRegressionNorm(X, Y, X, delta),'-',X,Y_true,'b');
grid
title(['Kernel Normal regression with delta = ' num2str(delta)])
xlabel ('x');ylabel('y, fit');
legend('y = cos(x) + x + eps', 'kernel reg.', 'y = cos(x)+x')
