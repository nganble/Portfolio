function [beta]=dis_Beta1(x,par)
global omega T
mysinc = @(x) sin(x)./x;

beta = mysinc((omega+x)*T).*(omega+x)*par.dt.*sin((omega+x)*par.dt)./(2*(1-cos((omega+x)*par.dt)))...
    + mysinc((omega-x)*T).*(omega-x)*par.dt.*sin((omega-x)*par.dt)./(2*(1-cos((omega-x)*par.dt)))...
    -(1/2)*mysinc(x*T).*(x*par.dt).*sin(x*par.dt)./(2*(1-cos(x*par.dt))); 

end