function [ cov_Non_Gauss ] = Psi( cov_Gauss, a, b, c, d )
cov_Non_Gauss=a*(b*cov_Gauss+c*cov_Gauss^2+d*cov_Gauss^3);
