function [ cov_Gauss ] = Gaussian_ACF( cov_exact, tau3i, tau4i, tau3j, tau4j )
gami=tau4i/0.1226;
h3i=9.21*tau3i./(11.68-2.5*gami);
h4i=(gami-1)./(11.68-2.5*gami);
ki=1./sqrt(1+2*h3i.^2+6*h4i+15*h4i.^2);
gamj=tau4j/0.1226;
h3j=9.21*tau3j./(11.68-2.5*gamj);
h4j=(gamj-1)./(11.68-2.5*gamj);
kj=1./sqrt(1+2*h3j.^2+6*h4j+15*h4j.^2);
%         if cov_exact(s,t)<-0.4698
%             cov_Gauss(s,t)=-1;
%             continue;
%         end
a=ki*kj;
b=1+3*h4i+3*h4j+9*h4i*h4j;
c=2*h3i*h3j;
d=6*h4i*h4j;
p=(3*b*d-c^2)/3/d^2;
q=2*c^3/(27*d^3)-b*c/(3*d^2)-cov_exact/(ki*kj*d);
delta=(p/3)^3+(q/2)^2;
r=sqrt(-p/3);
theta=acos(-q/2/r^3);
A=-q/2+sqrt(delta);B=-q/2-sqrt(delta);
if d<0
    cov_Gauss_min=max([-1 (-c+sqrt(c^2-3*b*d))/(3*d)]);
    cov_Gauss_max=min([1 (-c-sqrt(c^2-3*b*d))/(3*d)]);
    if cov_exact<Psi( cov_Gauss_min, a, b, c, d )
        cov_Gauss=cov_Gauss_min;
    end
    if cov_exact>Psi( cov_Gauss_max, a, b, c, d )
        cov_Gauss=cov_Gauss_max;
    end
    cov_Gauss=-2*r*cos((theta+pi)/3)-c/3/d;
else if delta>0
        cov_Gauss_min=-1;
        cov_Gauss_max=1;
        if cov_exact<Psi( cov_Gauss_min, a, b, c, d )
            cov_Gauss=cov_Gauss_min;
        end
        if cov_exact>Psi( cov_Gauss_max, a, b, c, d )
            cov_Gauss=cov_Gauss_max;
        end
        cov_Gauss=nthroot(A,3)+nthroot(B,3)-c/3/d;
else if c<0
        cov_Gauss_min=-1;
        cov_Gauss_max=min([1 (-c-sqrt(c^2-3*b*d))/(3*d)]);
        if cov_exact<Psi( cov_Gauss_min, a, b, c, d )
            cov_Gauss=cov_Gauss_min;
        end
        if cov_exact>Psi( cov_Gauss_max, a, b, c, d )
            cov_Gauss=cov_Gauss_max;
        end
        cov_Gauss=-2*r*cos((theta-pi)/3)-c/3/d;
else if c>0
        cov_Gauss_min=max([-1 (-c+sqrt(c^2-3*b*d))/(3*d)]);
        cov_Gauss_max=1;
        if cov_exact<Psi( cov_Gauss_min, a, b, c, d )
            cov_Gauss=cov_Gauss_min;
        end
        if cov_exact>Psi( cov_Gauss_max, a, b, c, d )
            cov_Gauss=cov_Gauss_max;
        end
        cov_Gauss=2*r*cos(theta/3)-c/3/d;
end
end
end
end
