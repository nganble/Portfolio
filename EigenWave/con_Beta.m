function [beta]=con_Beta(x)
global omega T
mysinc = @(u) sin(u)./u; % define sinc function
beta = mysinc((omega-x)*T)+mysinc((omega+x)*T)-(1/2)*mysinc(x*T);
end