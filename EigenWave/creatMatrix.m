function [A] = creatMatrix(par)
% Output: Matrix A of size Ngi*Ngi

% Simplifying notation
dr = par.dr; dth = par.dth;
i1r = par.i1r; i2r = par.i2r;
iath = par.iath; i2th = par.i2th;
i1th = par.i1th;
eqn = @(i,j) i+(par.Nr-1)*(j-1);
A = zeros(par.Ngi,par.Ngi);
for j = par.iath:par.i2th % j=0 to Nth
    for i = i1r:i2r % i=1 to Nr-1
        rij = par.r(i,j);
        k = eqn(i-1,j-1);
        A(k,k) = -2/dr^2-2/(rij*dth)^2; % Check
        if(j==iath) % Check
            A(k,eqn(i-1,i2th-1)) = 1/(rij*dth)^2;
            A(k,eqn(i-1,iath  )) = 1/(rij*dth)^2;
        elseif(j==i2th) % Check
            A(k,eqn(i-1,i1th-2)) = 1/(rij*dth)^2;
            A(k,eqn(i-1,i2th-2)) = 1/(rij*dth)^2;
        else % Check
            A(k,eqn(i-1,j-2)) = 1/(rij*dth)^2;
            A(k,eqn(i-1,j  )) = 1/(rij*dth)^2;
        end
        if(i==i1r) % Check
            A(k,eqn(i1r  ,j-1)) = 1/(dr^2)+1/(2*rij*dr);
        elseif(i==i2r) % Check
            A(k,eqn(i2r-2,j-1)) = 1/(dr^2)-1/(2*rij*dr);
        else % Check
            A(k,k-1) = 1/(dr^2)-1/(2*rij*dr);
            A(k,k+1) = 1/(dr^2)+1/(2*rij*dr);
        end
    end
end
end