function [tau3,tau4]=MCS_LM(x)
t=length(x);
x=sort(x);
b0=mean(x);
l(1)=b0;
b1=0;
for j=2:t
    b1=b1+(j-1)/(t-1)*x(j);
end
l(2)=2/t*b1-b0;
b2=0;
for j=3:t
    b2=b2+(j-1)*(j-2)/((t-1)*(t-2))*x(j);
end
l(3)=6/t*b2-6/t*b1+b0;
b3=0;
for j=4:t
    b3=b3+(j-1)*(j-2)*(j-3)/((t-1)*(t-2)*(t-3))*x(j);
end
l(4)=20/t*b3-30/t*b2+12/t*b1-b0;
tau3=l(3)/l(2);
tau4=l(4)/l(2);