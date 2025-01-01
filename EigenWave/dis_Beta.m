function [beta]=dis_Beta(x,par)
global omega T

beta=(1/2)*(cos(omega*0*par.dt)-1/4).*cos(x*0*par.dt);
for n=1:par.Nt-1
    beta=beta+(cos(omega*n*par.dt)-1/4).*cos(x*n*par.dt);
end
beta=beta+(1/2)*(cos(omega*par.Nt*par.dt)-1/4).*cos(x*par.Nt*par.dt);

beta=(2/T)*par.dt*beta;

end