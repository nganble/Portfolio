%------------------Eigenvalue problem: Lu+(lambda^2)u=0-----------------%
%------------------Dirichlet BCs----------------------------------------%
% Compute for 1D, 2D, and the annulus using both implicit and explicit scheme


clc;
close all;
clear all;
global omega T option scheme c numCount geometry cfl cpu;


geometry = '1D';
%geometry = '2D';
%geometry = 'Annulus';

dm = getDomain();

% Targeting eigenvalue
omega = 10;
Np = 1;
m = 2; % to adjust number of grid points
T = Np*(2*pi)/omega;
cfl = .6;
c = 1;
tFinal = T;
numCount = 0;


%option = 'testWaveSolver';
option = 'EigenWave';


scheme = 'Explicit';
%scheme = 'Implicit';

ms = 'eigen';
%ms = 'poly';
%ms = 'trig'; % Only for 2D and annulus
%ms = 'pulse'; % Only for 1D
%ms = 'mixed'; % Only for 1D

%methodName = "Power";
%methodName = "Arnoldi";
methodName = "eigs";

fun = getMS(dm,ms);

%%Test the waveSolver
if(strcmp(option,'testWaveSolver'))
    % Grid refinement
    numResolutions = 4;
    errMax = zeros(numResolutions,1);
    for m = 1:numResolutions
        par = setupMesh(dm,m);
        if (strcmp(geometry,'1D'))
            V0 = fun.u0(par.x);
            V1 = fun.u1(par.x);
            uexact = fun.ue(par.x,tFinal);
        elseif (strcmp(geometry,'2D')||strcmp(geometry,'Annulus'))
            if (strcmp(geometry,'Annulus')&&strcmp(ms,'eigen'))
                V0 = fun.u0(par.r,par.th);
                V1 = fun.u1(par.r,par.th);
                uexact = fun.ue(par.r,par.th,tFinal);
            else
                V0 = fun.u0(par.x,par.y);
                V1 = fun.u1(par.x,par.y);
                uexact = fun.ue(par.x,par.y,tFinal);
            end
        end
        if (strcmp(scheme,'Implicit'))
            [par.dA, par.decomTime] = formImplicitMatrix(par);
        end
        W = WaveHoltz(V0,V1,par,fun);
        errMax(m) = max(max(abs(uexact-W)));
        if (strcmp(geometry,'1D'))
            fprintf("%s:MS:%s: Scheme:%s m=%d: tFinal=%4.2f, Nx=%d\nTime to solve: %5.2e, maxErr=%7.2e",...
                geometry,ms,scheme,m,tFinal,par.Nx,cpu,errMax(m));
        elseif (strcmp(geometry,'2D'))
            fprintf("%s:MS:%s: Scheme:%s m=%d: tFinal=%4.2f Nx=Ny=%d Nt=%d\n",...
                geometry, ms, scheme, m, tFinal, par.Nx, par.Nt);
            if (strcmp(scheme, "Explicit"))
                fprintf("Solving time: %5.2f(s), maxErr=%7.2e", cpu,errMax(m));
            elseif (strcmp(scheme, "Implicit"))
                fprintf("decomTime=%5.2f(s) Solving time: %5.2f(s), maxErr=%7.2e", par.decomTime, cpu, errMax(m));
            end
        elseif (strcmp(geometry,'Annulus'))
            fprintf("%s:MS:%s: Scheme:%s m=%d: tFinal=%4.2f, Nr=%d, Nth=%d,\nTime to solve: %5.2e, maxErr=%7.2e",...
                geometry,ms,scheme,m,tFinal,par.Nr,par.Nth,cpu,errMax(m));
        end
        if(m==1)
            fprintf("\n");
        else
            fprintf(" ratio=%5.2f\n",errMax(m-1)/errMax(m));
        end
    end
    return;
end % End for test wave solver

%% EigenWave algorithm
% Set up grid points and grid functions
par = setupMesh(dm,m);
totalTime = 0;
% Initial guess
[W0,W1] = getICs(par);

% Discrete eigenvalues
if (strcmp(geometry,'1D'))
    lambda_h = zeros(par.Nx-1,1);
    for n=1:par.Nx-1
        lambda_h(n)=2*sin(n*pi*par.dx/2)/par.dx;
    end
elseif (strcmp(geometry,'2D'))
    lambda_h = zeros(par.Ngi,1);
    for iy = 1:par.Ny-1
        for ix = 1:par.Nx-1
            ie = par.eqn1(ix,iy);
            lambda_h(ie) = -4*(sin(ix*pi*par.dx/(2*dm.RB)))^2/(par.dx^2)...
                -4*(sin(iy*pi*par.dy/(2*dm.TB)))^2/(par.dy^2);
            lambda_h(ie) = sqrt(-lambda_h(ie));
        end
    end
elseif (strcmp(geometry,'Annulus'))
    L = creatMatrix(par);
    [VL,lambda_h] = eig(L,'vector');
%     % Check the residual
%     Res1 = zeros(par.Ngi,1);
%     for j = 1:par.Ngi
%         vj = VL(:,j);
%         vj = Vector_to_Matrix(vj,par);
%         Lhj = par.Lh(expand_matrix(vj,par),par.I,par.J);
%         vj = Matrix_to_Vector(vj);
%         Lhj = Matrix_to_Vector(Lhj);
%         Res1(j) = norm(Lhj-lambda_h(j)*vj,2);
%     end
%     resMax1 = max(abs(Res1));
%     fprintf("Checking residual of discrete matrix\nmaxRes = %5.2e\n",resMax1);
    lambda_h = sqrt(-lambda_h);
end

% Finding eigenvalues by differnt methods
if (strcmp(scheme,'Implicit'))
    [par.dA, par.decomTime] = formImplicitMatrix(par);
end
if(strcmp(methodName,"Power"))
    % Power method
    maxit = 2000;
    tol = 1e-8;
    
    W0 = W0/dis_norm(W0,par);
    
    numIter = 1; % to count total number of iterations
    % Compute convergence rate using eigenvector
    errOld = 1; % use this to initialize the ratio
    ratio1 = 1;
    % Compute convergence rate using residual
    resOld = 1;
    ratio2 = 1;
    comBeta = 0;
    
    %-----------------------Start power-iteration loop---------------------
    for h = 1:maxit
        % Apply WaveHoltz filter
        cpu0 = cputime;
        Z = WaveHoltz(W0,W1,par,fun); % Z=A*V
        totalTime = totalTime + cputime - cpu0;
        % Using Rayleigh quotient to find eigenvalue
        comBeta = dis_product(W0, Z, par);  % lambda=V'*A*V/||V||_2
        resNew = dis_norm(Z-comBeta*W0,par); % r=||A*V-lambda*V||_2
        ratio1 = resNew/resOld; % estimate the ratio using residue
        % flip sign if they are in the opposite direction
        if(comBeta<0)
            W0 = -W0;
        end
        Z = Z/dis_norm(Z,par); % normalize Z
        % Compute the 2-norm error
        errNew = dis_norm(Z-W0,par);
        ratio2 = errNew/errOld; % estimate the ratio using the errror of eigenvector
        % Print out the ratio after every 10 ierations
        if (mod(h,10)==0)
            fprintf("Iter %d: ratio_Res=%6.5f, Ratio_Vec=%6.5f\n",h,ratio1,ratio2);
        end
        % Stop the loop if the error exceeds the tolerance
        
        if (errNew<tol)
            break;
        end
        % Update for next iteration
        W0 = Z;
        errOld = errNew;
        resOld = resNew;
        numIter = numIter+1;
    end
    %-----------------------End power-iteration loop-----------------------

    comEigVec = W0;
    
%     % Eigenvalue for the original problem
%     Lhi = par.DpDm(W0,par.I);
%     Lhi = expand_matrix(Lhi,par);
%     comLambda = dis_product(W0,Lhi,par);
%     comLambda = sqrt(-comLambda);
    
elseif(strcmp(methodName,"Arnoldi"))
    % Arnoldi Method
    numTerm = 40;
    tol = 1e-7; % tolerance
    H = zeros(numTerm+1,numTerm);
    Q = zeros(par.NIt,(numTerm+1)*par.NJt);
    Q(:,1:par.NJt) = W0/dis_norm(W0,par);

    %-------------------Start loop for Arnoldi Algorithm------------------%
    for j = 1:numTerm
        cpu0 = cputime;
        v = WaveHoltz(Q(:,(j-1)*par.NJt+1:j*par.NJt),W1,par,fun);
        totalTime = totalTime + cputime - cpu0;
        for i=1:j
            H(i,j) = dis_product(Q(:,(i-1)*par.NJt+1:i*par.NJt),v,par);
            v = v-H(i,j)*Q(:,(i-1)*par.NJt+1:i*par.NJt);
        end
        H(j+1,j) = dis_norm(v,par);
        Q(:,j*par.NJt+1:(j+1)*par.NJt)=v/H(j+1,j);
    end
    
    %-------------------End loop for Arnoldi Algorithm--------------------%
    H = H(1:numTerm,:); % delete last row of H
    Q = Q(:,1:numTerm*par.NJt); % delete last column of Q
    % Convert matrix to vector form
    Q1 = zeros(par.Ng,numTerm);
    for j=1:numTerm
        Q1(:,j) = Matrix_to_Vector(Q(:,(j-1)*par.NJt+1:j*par.NJt));
    end
    [W,D] = eig(H,'vector');
    
    % Compute residue to find good approximation for eigpairs
    comBeta = zeros(numTerm,1);
    comEigVec = zeros(par.NIt,numTerm*par.NJt);
    count = 0;
    for j = 1:numTerm
        
        qj = Q1*W(:,j); % approximation of eigenvector
        beta = D(j); % approximation of eigenvalue
        qj = Vector_to_Matrix(qj,par); % convert qj to matrix form
        cpu0 = cputime;
        Lhj = WaveHoltz(qj,W1,par,fun);
        totalTime = totalTime + cputime - cpu0;
        res = dis_norm(Lhj-beta*qj,par)/abs(beta); %||A*q_j-lambda*q_j||_2
        if (res < tol)
            count = count+1;
            comBeta(count) = beta;
            %comEigVec(:,count) = qj;
            comEigVec(:,(count-1)*par.NJt+1:count*par.NJt) = qj;
        end
    end
    if(count==0)
        fprintf("\nNo eigenvalues found! Need to increase n.\n");
        return;
    else
        % Shrink the size
        comBeta = comBeta(1:count);
        comEigVec = comEigVec(:,1:count*par.NJt);
    end
    
    
elseif(strcmp(methodName,'eigs'))
    % Using built-in function eigs in MATLAB
    numRequest = 16;
    fe = @(x) myFunction(x,W1,par,fun);
    cpu0 = cputime;
    [comEigVec, comBeta] = eigs(fe,par.Ngi,numRequest);
    totalTime = totalTime + cputime - cpu0;
    comBeta = diag(comBeta);
end

% Compute eigenvalues of the original problem
numComputed = size(comBeta,1);
comLambda = zeros(numComputed,1);
for j=1:numComputed
    % Apply Rayleigh quotient lambda=V'*DpDm(V)/||V||_2
    if (strcmp(methodName,'eigs'))
        vj = Vector_to_Matrix(comEigVec(:,j),par);
        
        vj = expand_matrix(vj,par);
    else
        vj = comEigVec(:,(j-1)*par.NJt+1:j*par.NJt);
    end
    if (strcmp(geometry,'1D'))
        Lhj = par.Lh(vj,par.I);
    else
        Lhj = par.Lh(vj,par.I,par.J);
    end
    Lhj = expand_matrix(Lhj,par);
    comLambda(j) = dis_product(vj,Lhj,par)/dis_norm(vj,par)^2;
    comLambda(j) = sqrt(-comLambda(j));
end

% Find discrete eigenvalue corresponding to computed ones
disLambda = zeros(numComputed,1);
for j = 1:numComputed
    [val,index] = min(abs(comLambda(j)-lambda_h));
    disLambda(j) = lambda_h(index);
end

% Modified eigenvalues
if (strcmp(scheme,'Explicit'))
    modLambda = (2/par.dt)*asin(par.dt*disLambda/2);
elseif (strcmp(scheme,'Implicit'))
    modLambda = (2/par.dt)*asin(((1/2)*(disLambda*par.dt))./(sqrt(1+(1/2)*(disLambda*par.dt).^2)));
    
end
modBeta = dis_Beta(modLambda,par);


% Plot results
if (omega<21)
    numPoints=4*round(omega/pi)+1;
else
    numPoints=3*round(omega/pi)+1;
end
meshPoint = linspace(0,pi*numPoints,numPoints*50+1);
% Continuous WaveHolts time-filter
figure(1);
plot(meshPoint,con_Beta(meshPoint),'-b',"Linewidth",1.5);
hold on;
num = min(numPoints,par.Ngi);
plot(disLambda,con_Beta(disLambda),"*r","Linewidth",1.2);
hold on;
plot(comLambda,con_Beta(comLambda),"o k","Linewidth",1.2);
legend("\beta","discrte \lambda_j true","computed");
xlabel("\lambda");
ylabel("\beta");
title(sprintf("Filter function with omega=%d, Np=%d, %d eigs computed",...
    omega,Np,numComputed));
grid on;
% Discrete WaveHolts time-filter
figure(2);
plot(meshPoint,dis_Beta1(meshPoint,par),'-b','Linewidth',1.5);
hold on;
plot(modLambda,modBeta,'*r','Linewidth',1.2);
hold on;
plot(modLambda,comBeta,'o k','Linewidth',1.2);
legend("\beta","discrte \lambda_j","computed");
xlabel("\lambda");
ylabel("\beta");
title(sprintf("Discrete Filter function, %d eigs computed",numComputed));
grid on;

% Print out results
comLambda = sort(comLambda,'descend');
disLambda = sort(disLambda,'descend');
errMax = abs(comLambda-disLambda);
if (strcmp(geometry, '1D'))
    fprintf("EigenWave algorithm: Nx=%d Nt=%d tFinal=%4.2f dx=%4.2e dt=%4.2e\n",...
        par.Nx, par.Nt, T, par.dx, par.dt);
elseif (strcmp(geometry, '2D'))

    fprintf("EigenWave algorithm: Nx=Ny=%d Nt=%d tFinal=%4.2f dx=%4.2e dy=%4.2e dt=%4.2e\n",...
        par.Nx, par.Nt, T, par.dx, par.dy, par.dt);
else
    fprintf("EigenWave algorithm: Nr=%d Nth=%d Nt=%d tFinal=%4.2f dr=%4.2e dth=%4.2e dt=%4.2e\n",...
        par.Nr, par.Nth, par.Nt, T, par.dr, par.dth, par.dt);

end
if (strcmp(scheme,'Implicit'))
    fprintf('decomTime:%5.2f(s) totalTime=%5.2f(s)\n', par.decomTime, totalTime);
elseif (strcmp(scheme,'Explicit'))
    fprintf('totalTime=%5.2f(s)\n', totalTime);
end
fprintf('Domain:%s, Time_per_iteration: %5.2e\n',geometry, totalTime/numCount);
fprintf('Scheme:%s, Method=%s\nnumComputed = %d,  numWaveHoltzSolver = %d\n',...
    scheme,methodName,numComputed,numCount);
fprintf("    disLambda           comLambda          err\n");
for j = 1:numComputed
    fprintf("%15.10f     %15.10f      %5.2e\n",...
        disLambda(j),comLambda(j),errMax(j));
end
