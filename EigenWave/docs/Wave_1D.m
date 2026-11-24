%---------------Solving w_tt=c^2w_xx+f(x,t),           a<x<b, 0<t<tFinal
%---------------------- w(a,t)=ga(t),   w(b,t)=gb(t),  0<=t<=tFinal 
%---------------------- w(x,0)=v(x),  w_t(x,0)=v1(x),  a<=x<=b
% Explicit and implicit schemes using second-order finite difference
% approximations (see waveSchemes.pdf, eq. (2) and eq. (8)-(18))

clc;
close all;
clear;

% Innitial set-up
a=0; b=1; c=1;
cfl=.9;
tFinal=.65;

scheme='Explicit';
% scheme='Implicit';

% ms='poly';
% ms='pulse';
ms='eigen';
if(strcmp(ms,'poly'))
    % Polynomial manufactured solution
    b0=.5; b1=.7; b2=.9; % time coeffivelcients
    c0=1; c1=.8; c2=.6;  % space coefficients
    ue=@(x,t) (b0+b1*t+b2*t^2)*(c0+c1*x+c2*x.^2);
    ut=@(x,t) (b1+2*b2*t)*(c0+c1*x+c2*x.^2);
    utt=@(x,t) (2*b2)*(c0+c1*x+c2*x.^2);
    uxx=@(x,t) (b0+b1*t+b2*t^2).*2*c2;
elseif(strcmp(ms,'pulse'))
    % Traveling pulse
    beta=20; x0=.25;
    ue=@(x,t) exp(-(beta^2)*(x-x0-c*t).^2);
    ut=@(x,t) 2*c*beta^2*(x-x0-c*t).*ue(x,t);
    utt=@(x,t) (4*c^2*beta^4*(x-x0-c*t).^2-2*c^2*beta^2).*ue(x,t);
    ux=@(x,t) -2*beta^2*(x-x0-c*t).*ue(x,t);
    uxx=@(x,t) (4*beta^4*(x-x0-c*t).^2-2*beta^2).*ue(x,t);
else
    % Test Eigenfunction
    n=1;
    lambda=n*pi;
    ue=@(x,t) cos(lambda*c*t)*sin(lambda*x);
    ut=@(x,t) -lambda*c*sin(lambda*c*t)*sin(lambda*x);
    utt=@(x,t) -lambda^2*c^2*cos(lambda*c*t)*sin(lambda*x);
    ux=@(x,t) lambda*cos(lambda*c*t)*cos(lambda*x);
    uxx=@(x,t) -lambda^2*cos(lambda*c*t)*sin(lambda*x);
    
end
% Forcing term
f=@(x,t) utt(x,t)-c^2*uxx(x,t);
% Boundary condition
ga=@(t) ue(a,t);
gb=@(t) ue(b,t);
% Initial condition
u0=@(x) ue(x,0);
u1=@(x) ut(x,0);
u0xx=@(x) uxx(x,0);
numResolutions=6;
errMax=zeros(numResolutions,1);
errL1=zeros(numResolutions,1);

for m=1:numResolutions
    Nx=20*2^m;
    dx=(b-a)/Nx;
    x=linspace(a,b,Nx+1)';
    ia=1; %index of boundary point at x=a
    ib=ia+Nx; %index of boundary point at x=b
    i1=ia+1; %first interior point
    i2=ib-1; %last interior point
    I=ia:ib;
    I1=i1:i2;
   
    
    % D+D- operator
    DpDm=@(u,I) (u(I+1)-2*u(I)+u(I-1))/(dx^2);
    
    % Time-stepping
    dt=cfl*dx/c;
    Nt=ceil(tFinal/dt);
    dt=tFinal/Nt;

    % Implicit matrix, eq. (11)
    if(strcmp(scheme,'Implicit'))
        r=(c^2*dt^2)/(2*dx^2);
        A=spalloc(Nx+1,Nx+1,3*(Nx+1));
        A(ia,ia)=1; % BC at x = a
        for i=i1:i2
            A(i,i-1)=-r;
            A(i,i  )=1+2*r;
            A(i,i+1)=-r;
        end
        A(ib,ib)=1; % BC at x = b
        dA=decomposition(A); % factor once, reuse every step
    end

    % Allocate space for solution at three levels
    unm1=zeros(Nx+1,1);
    un=zeros(Nx+1,1);
    unp1=zeros(Nx+1,1);
    
    t=0;
    unm1=u0(x);
    % Old first step: used the exact second derivative u0xx and no BCs
    % un=u0(x)+dt*u1(x)+(dt^2/2)*(c^2*u0xx(x)+f(x,0));
    if(strcmp(scheme,'Explicit'))
        % First step, eq. (2c): Taylor expansion with discrete D+D- of V0
        un=u0(x)+dt*u1(x);
        un(I1)=un(I1)+(dt^2/2)*(c^2*DpDm(unm1,I1)+f(x(I1),0));
        un(ia)=ga(dt); % BC at x = a
        un(ib)=gb(dt); % BC at x = b
    elseif(strcmp(scheme,'Implicit'))
        % First step, eq. (17)-(18): eliminate W^{-1} using D0t W^0 = V1
        V0=u0(x);
        V1=u1(x);
        rhs=zeros(Nx+1,1);
        rhs(I1)=V0(I1)+dt*V1(I1)-(c^2*dt^3/2)*DpDm(V1,I1)...
            +(dt^2/4)*(f(x(I1),dt)+f(x(I1),-dt));
        rhs(ia)=ga(dt); % BC at x = a
        rhs(ib)=gb(dt); % BC at x = b
        un=dA\rhs;
    end

    %--------Start time-stepping loop--------%
    for n=1:Nt-1

        t=n*dt; % new time
        if(strcmp(scheme,'Explicit'))
            % Eq. (2a)
            unp1(I1)=2*un(I1)-unm1(I1)+(c^2*dt^2)*DpDm(un,I1)+(dt^2)*f(x(I1),t);
            % Update boundary point
            unp1(ia)=ga(t+dt);
            unp1(ib)=gb(t+dt);
        elseif(strcmp(scheme,'Implicit'))
            % Eq. (10)-(12): A*W^{n+1} = G^n
            rhs=zeros(Nx+1,1);
            rhs(I1)=2*un(I1)-unm1(I1)+(c^2*dt^2/2)*DpDm(unm1,I1)...
                +(dt^2/2)*(f(x(I1),t+dt)+f(x(I1),t-dt));
            rhs(ia)=ga(t+dt); % BC at x = a
            rhs(ib)=gb(t+dt); % BC at x = b
            unp1=dA\rhs;
        end

%         % Make movie
%         plot(x,unp1);
%         grid on;
%         drawnow;
        
        % Update unm1 and un for next step
        unm1=un;
        un=unp1;  
    end
    %%%%%-----------End time-stepping loop----------%%%%%
    
    % Compute errors
    t=tFinal; % tfinal time
    uexact=ue(x,t); % evaluate exact solution at final time
    err=uexact-un;
    errMax(m)=max(abs(err));
    fprintf('ms=%s: t=%8.2e: m=%d Nx=%4d Nt=%4d dt=%8.2e maxErr=%8.2e',...
        ms,t,m,Nx,Nt,dt,errMax(m));
    if(m==1)
        fprintf('\n');
    else
        fprintf(' ratio=%5.2f rate=%5.2f\n',errMax(m-1)/errMax(m),log2(errMax(m-1)/errMax(m)));
    end
    
    
    % PLot results
    if((strcmp(ms,'pulse'))||(strcmp(ms,'eigen')))
        if(m<=4)
            figure();
            plot(x,un,'b-o',x,uexact,'k','Linewidth',1.5);
            legend('computed', 'true','location','northeast');
            title(sprintf('Wave Equation:ms:%s, Nx=%d, t=%2.2f',...
                ms,Nx,t));
            xlabel('x'); ylabel('u');
            grid on;
            % Plot error
            figure();
            plot(x,err,'r-+','Linewidth',1.5);
            legend('Error');
            title(sprintf('Error:ms:%s, Nx=%d, t=%2.2f',...
                ms,Nx,t));
            grid on;
        end
    end
end
