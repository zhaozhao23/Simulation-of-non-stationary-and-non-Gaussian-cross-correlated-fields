function [ tao3, tao4 ] = L_moments_stat( )
%This function is used to compute the target L-skewness/kurtosis of Beta distribution
f2=@(F) (betainv(F,4,2)*(1.87+3.74)-3.74).*(2*F-1);
lamda2=integral(f2,0,1);
f3=@(F) (betainv(F,4,2)*(1.87+3.74)-3.74).*(6*F.^2-6*F+1);
lamda3=integral(f3,0,1);
f4=@(F) (betainv(F,4,2)*(1.87+3.74)-3.74).*(20*F.^3-30*F.^2+12*F-1);
lamda4=integral(f4,0,1);
tao3=lamda3/lamda2;
tao4=lamda4/lamda2;    
    


