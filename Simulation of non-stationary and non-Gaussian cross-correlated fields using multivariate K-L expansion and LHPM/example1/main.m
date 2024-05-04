clear all;clc;close all;
tic;
a11=1/10;
a22=1/5;
a12=(a11+a22)/2;
C11=@(t,s) exp(-a11*abs(t-s));
C22=@(t,s) exp(-a22*abs(t-s));
C12=@(t,s) 0.7*exp(-a12*abs(t-s));
tau3=-0.0989;tau4=0.0814;%Target L-skewness/kurtosis
L=10;q=50;
dx=L/q;
x=[dx/2:dx:L-dx/2];

Sig=zeros(2*q,2*q);
for i=1:q
    for j=1:q
        Sig((i-1)*2+1,(j-1)*2+1)=Gaussian_ACF( C11(x(i),x(j)), tau3, tau4, tau3, tau4 );
        Sig((i-1)*2+1,(j-1)*2+2)=Gaussian_ACF( C12(x(i),x(j)), tau3, tau4, tau3, tau4 );
        Sig((i-1)*2+2,(j-1)*2+1)=Gaussian_ACF( C12(x(i),x(j)), tau3, tau4, tau3, tau4 );
        Sig((i-1)*2+2,(j-1)*2+2)=Gaussian_ACF( C22(x(i),x(j)), tau3, tau4, tau3, tau4 );
    end
end
[V,D]=eig(Sig);
[D_sort,index]=sort(diag(D),'descend');
lamda=D_sort*dx;
V_sort=V(:,index);
M=12;
dt=0.1;
t=[0:dt:10];
nt=length(t);
Nsim=1e5;
xi=randn(Nsim,M);
U=zeros(Nsim,2,nt);
Y=zeros(Nsim,2,nt);
for i=1:nt
    for k=1:q
        Sigqi((k-1)*2+1,1)=Gaussian_ACF( C11(t(i),x(k)), tau3, tau4, tau3, tau4 );
        Sigqi((k-1)*2+1,2)=Gaussian_ACF( C12(t(i),x(k)), tau3, tau4, tau3, tau4 );
        Sigqi((k-1)*2+2,1)=Gaussian_ACF( C12(t(i),x(k)), tau3, tau4, tau3, tau4 );
        Sigqi((k-1)*2+2,2)=Gaussian_ACF( C22(t(i),x(k)), tau3, tau4, tau3, tau4 );
    end
    for j=1:M
        eigfun(j,i,:)=1/lamda(j)*sqrt(dx)*Sigqi'*V_sort(:,j);
        U(:,:,i)=U(:,:,i)+sqrt(lamda(j))*xi(:,j)*squeeze(eigfun(j,i,:))';
    end
    Y(:,1,i)=LHPM(U(:,1,i),tau3,tau4);
    Y(:,2,i)=LHPM(U(:,2,i),tau3,tau4);
end
toc;
%Estimated L-skewness/kurtosis from the generated non-stationary non-Gaussian samples
for i=1:nt
    [tau3y1s(i),tau4y1s(i)]=MCS_LM(Y(:,1,i));
    [tau3y2s(i),tau4y2s(i)]=MCS_LM(Y(:,2,i));
end

figure;
plot(t,tau3*ones(1,nt),'r--','LineWidth', 1);hold on;
plot(t,tau3y1s,'k');
plot(t,tau3y2s,'b-.');
legend('Target','$Y_1(x)$', '$Y_2(x)$', 'Interpreter', 'latex');
xlabel('\itx\rm');
ylabel('L-Skewness');
set(gca,'FontSize',14);
set(gca,'FontName','Times New Roman');
ylim([-0.4 0.2]);
figure;
plot(t,tau4*ones(1,nt),'r--','LineWidth', 1);hold on;
plot(t,tau4y1s,'k');
plot(t,tau4y2s,'b-.');
legend('Target','$Y_1(x)$', '$Y_2(x)$', 'Interpreter', 'latex');
xlabel('\itx\rm');
ylabel('L-Kurtosis');
set(gca,'FontSize',14);
set(gca,'FontName','Times New Roman');
ylim([-0.2 0.4]);

figure;
plot(t,squeeze(Y(1:100,1,:)));
xlabel('$x$', 'Interpreter', 'latex');ylabel('$Y_1(x)$', 'Interpreter', 'latex');
set(gcf,'color','w');
figure;
plot(t,squeeze(Y(1:100,2,:)));
xlabel('$x$', 'Interpreter', 'latex');ylabel('$Y_2(x)$', 'Interpreter', 'latex');
set(gcf,'color','w');
%Target covariance function
for i=1:nt
    for j=1:nt
        rohy1y1t(i,j)=C11(t(i),t(j));
        rohy1y2t(i,j)=C12(t(i),t(j));
        rohy2y2t(i,j)=C22(t(i),t(j));
    end
end
%Estimated covariance function from the generated non-stationary non-Gaussian samples
for i=1:nt
    for j=1:nt
        rohy1y1s(i,j)=mean(Y(:,1,i).*Y(:,1,j));
        rohy1y2s(i,j)=mean(Y(:,1,i).*Y(:,2,j));
        rohy2y2s(i,j)=mean(Y(:,2,i).*Y(:,2,j));
    end
end

figure;
mesh(t,t,rohy1y1t); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('$C_{Y_1Y_1}(x,x^{\prime})$', 'Interpreter', 'latex');
set(gcf,'color','w');
figure;
mesh(t,t,rohy1y1s); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('$C_{Y_1Y_1}(x,x^{\prime})$', 'Interpreter', 'latex');
set(gcf,'color','w');zlim([min(min(rohy1y1s)) 1])
figure;
mesh(t,t,abs(1-rohy1y1s./rohy1y1t)); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('Error');
set(gcf,'color','w');
figure;
mesh(t,t,rohy1y2t); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('$C_{Y_1Y_2}(x,x^{\prime})$', 'Interpreter', 'latex');
set(gcf,'color','w');
figure;
mesh(t,t,rohy1y2s); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('$C_{Y_1Y_2}(x,x^{\prime})$', 'Interpreter', 'latex');
set(gcf,'color','w');
figure;
mesh(t,t,abs(1-rohy1y2s./rohy1y2t)); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('Error');
set(gcf,'color','w');
figure;
mesh(t,t,rohy2y2t); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('$C_{Y_2Y_2}(x,x^{\prime})$', 'Interpreter', 'latex');
set(gcf,'color','w');
figure;
mesh(t,t,rohy2y2s); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('$C_{Y_2Y_2}(x,x^{\prime})$', 'Interpreter', 'latex');
set(gcf,'color','w');zlim([min(min(rohy2y2s)) 1])
figure;
mesh(t,t,abs(1-rohy2y2s./rohy2y2t)); colormap bone; box on; axis tight;
xlabel('$x$', 'Interpreter', 'latex');ylabel('$x^{\prime}$', 'Interpreter', 'latex');zlabel('Error');
set(gcf,'color','w');

eigfunx1=squeeze(eigfun(:,:,1));
figure;
hold on;
for i=1:12
    plot(t,eigfunx1(i,:),'LineWidth',1);
end
hold off;
xlabel('$x$','Interpreter','latex');
ylabel('$\phi_{1i}$','Interpreter','latex');
set(gca,'FontSize',14);
set(gca,'FontName','Times New Roman');

eigfunx2=squeeze(eigfun(:,:,2));
figure;
hold on;
for i=1:12
    plot(t,eigfunx2(i,:),'LineWidth',1);
end
hold off;
xlabel('$x$','Interpreter','latex');
ylabel('$\phi_{2i}$','Interpreter','latex');
set(gca,'FontSize',14);
set(gca,'FontName','Times New Roman');

figure;
yy=-3:0.1:2;
ff=betacdf((yy+3.74)/(1.87+3.74),4,2);
plot(yy,ff,'r--','LineWidth', 1);
hold on;
[aa,bb]=ecdf(Y(:,1,i));
plot(bb,aa,'b');
legend('Target','Estimated');
xlabel('$Y_1(x)$', 'Interpreter', 'latex');
ylabel('CDF');
set(gca,'FontSize',14);
set(gca,'FontName','Times New Roman');
figure;
yy=-3:0.1:2;
ff=betacdf((yy+3.74)/(1.87+3.74),4,2);
plot(yy,ff,'r--','LineWidth', 1);
hold on;
[aa,bb]=ecdf(Y(:,2,i));
plot(bb,aa,'b');
legend('Target','Estimated');
xlabel('$Y_2(x)$', 'Interpreter', 'latex');
ylabel('CDF');
set(gca,'FontSize',14);
set(gca,'FontName','Times New Roman');