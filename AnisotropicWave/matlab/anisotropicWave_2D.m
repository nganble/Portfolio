%--------------------- Anisotropic wave equationin 2D --------------------%
%--------------------- rho*u_tt=\div(m*\div u) ---------------------------%
%--------------------- u=u(x,y,t), rho=rho(x,y)---------------------------%
%---------------------          [ m_11(x,y) m_12(x,y)]--------------------%
%--------------------- m=m(x,y)=[                    ]--------------------%
%---------------------          [ m_21(x,y) m_22(x,y)]--------------------%
%--------------------- (x,y)\in [a,b]*[c,d], t>0 -------------------------%




if ~exist('cases','var')        % standalone run: set the case here
    clc;
    ms      = "gauss";           % 'poly' | 'trig' | 'pulse' | 'trigPulse' | 'gauss' | 'gaussForce'
    coeffs  = "piecewise";          % 'const' | 'var' | 'piecewise'
    saveRef = false;            % true: write reference, false: compare
end
clearvars -except ms coeffs saveRef cases nCases k pass
close all;

%%  Domain
ax = -1; bx = 1;
ay = -1; by = 1;

% Wavenumbers
kx = 1;
ky = 1;


plotSol = false;
makeMovie = 0;  % 0=no movie (convergence test) | 1=make movie | 2=save movie too
topView = "on"; % 'on' | 'off'
if (makeMovie==0)
    tFinal=.5;
else
    if (strcmp(coeffs,"const"))
        if (strcmp(ms, "trig"))
            tFinal = 25;
        elseif (strcmp(ms, "pulse"))
            tFinal = 4;
        elseif (strcmp(ms, "trigPulse"))
            tFinal = 5;
        elseif (strcmp(ms, "gaussForce"))
            tFinal = 20;
        elseif (strcmp(ms, "gauss"))
            tFinal = 10;
        else
            tFinal = 10;
        end
    else
        tFinal = 1.2;
    end
    
end

cfl = 0.9;
debug = 0; % debug=1 check after first time-step | debug=2 check approx U1

% Grid refinement
if (makeMovie==0)
    startInd = 1;
    numResolutions = 4;
else
    if (strcmp(coeffs,"const"))
        if (strcmp(ms,"gauss"))
            startInd = 3;
            numResolutions = 3;
        elseif (strcmp(ms,"gaussForce"))
            startInd = 3;
            numResolutions = 3;
        else
            startInd = 3;
            numResolutions = 3;
        end
    else
        if (strcmp(ms,"gauss"))
            startInd = 3;
            numResolutions = 3;
        elseif (strcmp(ms,"gaussForce"))
            startInd = 4;
            numResolutions = 4;
        else
            startInd = 4;
            numResolutions = 4;
        end
    end
end


%%  Function rho and matrix m
if (strcmp(coeffs,"const"))            % isotropic: rho = 1, m = I
    RHO = @(x,y) 1 + 0*x + 0*y;
    M11 = @(x,y) 1 + 0*x + 0*y;   M11x = @(x,y) 0*x + 0*y;
    M12 = @(x,y) 0*x + 0*y;       M12x = @(x,y) 0*x + 0*y;
    M21 = @(x,y) 0*x + 0*y;       M21y = @(x,y) 0*x + 0*y;
    M22 = @(x,y) 1 + 0*x + 0*y;   M22y = @(x,y) 0*x + 0*y;
elseif (strcmp(coeffs,"var"))          % symmetric positive definite on (-1,1)^2
    RHO = @(x,y) 0.1*(sin(x+y)+1.1);
    M11 = @(x,y) 0.3*(x+y) + 2;     M11x = @(x,y) 0.3 + 0*x + 0*y;
    M12 = @(x,y) 0.5 + 0.2*(x+y);   M12x = @(x,y) 0.2 + 0*x + 0*y;
    M21 = M12;                      M21y = @(x,y) 0.2 + 0*x + 0*y;
    M22 = @(x,y) 0.6*(2*x+y) + 3;   M22y = @(x,y) 0.6 + 0*x + 0*y;

    % Another option for variable density and stiffness matrix
    % RHO = @(x,y) 0.1*(sin(x+y)+1.1);
    % M11 = @(x,y) 2 + 0.3*(x+y);                         M11x = @(x,y) 0.3 + 0*x + 0*y;
    % M12 = @(x,y) 0.5 + 0.2*sin(pi*x).*cos(pi*y);        M12x = @(x,y)  0.2*pi*cos(pi*x).*cos(pi*y);
    % M21 = M12;                                          M21y = @(x,y) -0.2*pi*sin(pi*x).*sin(pi*y);
    % M22 = @(x,y) 3 + 0.6*(2*x+y);                       M22y = @(x,y) 0.6 + 0*x + 0*y;
elseif (strcmp(coeffs,"piecewise"))    % left: x<0, right: x>=0
    L = @(x) (x < 0);  R = @(x) (x >= 0);
    RHO = @(x,y) (1*L(x) + 2*R(x)) + 0*y;
    M11 = @(x,y) (1*L(x) + 4*R(x)) + 0*y;
    M12 = @(x,y) (0*L(x) + 0.5*R(x)) + 0*y;
    M21 = @(x,y) (0*L(x) + 0.5*R(x)) + 0*y;
    M22 = @(x,y) (1*L(x) + 2*R(x)) + 0*y;
    % derivatives are not used: piecewise runs only with "gauss" (f = 0)
end

% Dispersive relation -- only use for constant variable cases
omega = sqrt(M11(ax,ay)*kx^2+(M12(ax,ay)+M21(ax,ay))*kx*ky+...
        M22(ax,ay)*ky^2)/sqrt(RHO(ax,ay));

%% Manufactured and exact solution for testing
if (strcmp(ms,'poly'))
    % Polynomial manufactured solution
    a0 = .5; a1 = .7; a2 = .9;     % time coefficients
    b0 = .5; b1 = .8; b2 = .6;     % space coefficients in x
    c0 = .3; c1 = .25; c2 = .2;    % space coefficients in y

    ue = @(x,y,t) (b0+x.*(b1+x*b2)).*(c0+y.*(c1+y*c2))*(a0+t*(a1+t*a2));
    uet = @(x,y,t) (b0+x.*(b1+x*b2)).*(c0+y.*(c1+y*c2))*(a1+2.*t*a2);
    uex = @(x,y,t) (b1+2.*x*b2).*(c0+y.*(c1+y*c2))*(a0+t*(a1+t*a2));
    uey = @(x,y,t) (b0+x.*(b1+x*b2)).*(c1+2.*y*c2)*(a0+t*(a1+t*a2));
    uexy = @(x,y,t) (b1+2.*x*b2).*(c1+2.*y*c2)*(a0+t*(a1+t*a2));
    uett = @(x,y,t) (b0+x.*(b1+x*b2)).*(c0+y.*(c1+y*c2))*(2.*a2+0.*t);
    uexx = @(x,y,t) (2.*b2+0.*x).*(c0+y.*(c1+y*c2))*(a0+t*(a1+t*a2));
    ueyy = @(x,y,t) (b0+x.*(b1+x*b2)).*(2.*c2+0.*y)*(a0+t*(a1+t*a2));
elseif (strcmp(ms,"trig"))
    % Exact solution
    ue  = @(x,y,t) 3*cos(kx*x+ky*y-omega*t);
    uet = @(x,y,t) 3*omega*sin(kx*x+ky*y-omega*t);
    %ue = @(x,y,t) 3*sin(kx*x+ky*y-omega*t)+3*cos(kx*x+ky*y-omega*t);
    %uet = @(x,y,t) -3*omega*cos(kx*x+ky*y-omega*t)+3*omega*sin(kx*x+ky*y-omega*t);
elseif (strcmp(ms,"pulse"))
    beta = 3;
    x0 = 0; % Central of the pulse
    ue = @(x,y,t) 4*exp(-beta*(kx*x+ky*y-omega*t-x0).^2);
    uet = @(x,y,t) 2*beta*omega*(kx*x+ky*y-omega*t-x0).*ue(x,y,t);
elseif (strcmp(ms,"trigPulse"))
    beta = 2.5;
    x0 = 0;
    kn = 2;
    ue = @(x,y,t) 2*exp(-beta*(kx*x+ky*y-omega*t-x0).^2)...
            .*cos(kn*(kx*x+ky*y-omega*t-x0));
    uet = @(x,y,t) 2*beta*omega*(kx*x+ky*y-omega*t-x0).*ue(x,y,t)...
            +kn*omega*2*exp(-beta*(kx*x+ky*y-omega*t-x0).^2)...
            .*sin(kn*(kx*x+ky*y-omega*t-x0));
    
end

% Initial conditions and boundary conditions and forcing
if (strcmp(ms,"gauss"))
    x0 = -0.5; y0 = (ay+by)/2;
    gamma = 1; sigmaX = 0.15; sigmaY = 0.15;
    u0 = @(x,y) exp(-gamma*(((x-x0)/sigmaX).^2+((y-y0)/sigmaY).^2));
    u1 = @(x,y) 0.*x+0.*y;
    gax = @(y,t) 0.*t+0.*y;
    gbx = @(y,t) 0.*t+0.*y;
    gay = @(x,t) 0.*x+0.*t;
    gby = @(x,t) 0.*x+0.*t;
    f = @(x,y,t) 0.*x+0.*y+0*t;
elseif (strcmp(ms,"gaussForce"))
    x0 = -0.5; y0 = (ay+by)/2;
    gamma = 3; sigmaX = 0.15; sigmaY = 0.15;
    u0 = @(x,y) 0.*x+0.*y;
    u1 = @(x,y) 0.*x+0.*y;
    gax = @(y,t) 0.*t+0.*y;
    gbx = @(y,t) 0.*t+0.*y;
    gay = @(x,t) 0.*x+0.*t;
    gby = @(x,t) 0.*x+0.*t;
    f = @(x,y,t) 2*exp(-gamma*(((x-x0)/sigmaX).^2+((y-y0)/sigmaY).^2)).*sin(gamma*t);
else
    
    u0   = @(x,y) ue(x,y,0);
    u1   = @(x,y) uet(x,y,0);
    u0xx = @(x,y) uedxdx(x,y,0);
    u0yy = @(x,y) uedydy(x,y,0);
    u1xx = @(x,y) uetxx(x,y,0);
    u1yy = @(x,y) uetyy(x,y,0);
    
    % Boundary conditions
    gax = @(y,t) ue(ax,y,t);
    gbx = @(y,t) ue(bx,y,t);
    gay = @(x,t) ue(x,ay,t);
    gby = @(x,t) ue(x,by,t);
    
    % Forcing term
    if (strcmp(ms, "poly"))
        f = @(x,y,t) RHO(x,y).*uett(x,y,t)-M11(x,y).*uexx(x,y,t)...
            -M22(x,y).*ueyy(x,y,t)-(M12(x,y)+M21(x,y)).*uexy(x,y,t)...
            -(M11x(x,y)+M21y(x,y)).*uex(x,y,t)...
            -(M12x(x,y)+M22y(x,y)).*uey(x,y,t);
    else
        f = @(x,y,t) 0.*x+0.*y+0*t;
    end
end

if (makeMovie>1)
    % Define the parameters of the video
    videoFile = 'myVideo.mp4';
    fps = 40; % frames per second
    % Create a VideoWriter object
    videoWriterObj = VideoWriter(videoFile, 'MPEG-4');
    videoWriterObj.FrameRate = fps;
    % Open the video file for writing
    open(videoWriterObj);
end

%% Convergence test
errMax = zeros(numResolutions,1);
NxList = zeros(numResolutions,1); NtList = NxList; dtList = NxList; relErr = NxList;
for m = startInd:numResolutions
    Nx = 10*2^m; Ny = 10*2^m;
    dx = (bx-ax)/Nx;
    dy = (by-ay)/Ny;
    iax = 1; iay = 1; 
    ibx = iax+Nx; iby = iay+Ny;
    i1x = iax+1;  i1y = iay+1;
    i2x = ibx-1;  i2y = iby-1;
    I  = i1x:i2x; % interior points in x
    I1b = iax:ibx; % interior points & boundary points in x
    J  = i1y:i2y; % interior points in y
    I2b = iay:iby; % interior points & boundary points in y
    Ngx = Nx+1; % total grid points in x
    Ngy = Ny+1; % total grid points in y
    x = zeros(Nx+1,Ny+1);
    y = zeros(Nx+1,Ny+1);
    % fill in x&y
    for iy = 1:Ngy
        for ix = 1:Ngx
            x(ix,iy) = ax+(ix-iax)*dx;
            y(ix,iy) = ay+(iy-iay)*dy;
        end
    end
    % Matrices hold grid values for rho and m
    rho = RHO(x,y);
    m11 = M11(x,y);
    m12 = M12(x,y);
    m21 = M21(x,y);
    m22 = M22(x,y);
    
%     % Define difference operators
%     xDef = @(u,m,I,J) (m(x(I,J)+dx/2,y(I,J)).*u(I+1,J)...
%         -(m(x(I,J)+dx/2,y(I,J))+m(x(I,J)-dx/2,y(I,J))).*u(I,J)...
%         +m(x(I,J)-dx/2,y(I,J)).*u(I-1,J))/(dx^2);
%     yDef = @(u,m,I,J) (m(x(I,J),y(I,J)+dy/2).*u(I,J+1)...
%         -(m(x(I,J),y(I,J)+dy/2)+m(x(I,J),y(I,J)-dy/2)).*u(I,J)...
%         +m(x(I,J),y(I,J)-dy/2).*u(I,J-1))/(dy^2);
%     xyDef = @(u,m,I,J) (m(x(I+1,J),y(I,J)).*(u(I+1,J+1)-u(I+1,J-1))...
%         -m(x(I-1,J),y(I,J)).*(u(I-1,J+1)-u(I-1,J-1)))/(4*dx*dy);
%     yxDef = @(u,m,I,J) (m(x(I,J),y(I,J+1)).*(u(I+1,J+1)-u(I-1,J+1))...
%         -m(x(I,J),y(I,J-1)).*(u(I+1,J-1)-u(I-1,J-1)))/(4*dx*dy);
    
    % New difference operators
    xDef = @(u,m,I,J) ((1/2)*(m(I+1,J)+m(I,J)).*u(I+1,J)...
        -((1/2)*(m(I+1,J)+2*m(I,J)+m(I-1,J))).*u(I,J)...
        +((1/2)*(m(I,J)+m(I-1,J))).*u(I-1,J))/(dx^2);
    yDef = @(u,m,I,J) ((1/2)*(m(I,J+1)+m(I,J)).*u(I,J+1)...
        -((1/2)*(m(I,J+1)+2*m(I,J)+m(I,J-1))).*u(I,J)...
        +((1/2)*(m(I,J)+m(I,J-1))).*u(I,J-1))/(dy^2);
    xyDef = @(u,m,I,J) ( m(I+1,J).*(u(I+1,J+1)-u(I+1,J-1))...
        -m(I-1,J).*(u(I-1,J+1)-u(I-1,J-1)) )/(4*dx*dy);
    yxDef = @(u,m,I,J) (m(I,J+1).*(u(I+1,J+1)-u(I-1,J+1))...
        -m(I,J-1).*(u(I+1,J-1)-u(I-1,J-1)))/(4*dx*dy);
    
    % Time stepping
    dt = cfl*min(min(sqrt(rho)./sqrt( m11/(dx^2)+abs((m12...
        +m21))/(4*dx*dy)+m22/(dy^2) )));
    Nt = ceil(tFinal/dt);
    dt = tFinal/Nt; % adjust dt to exactly reach tFinal
    
    % Allocate space for solution
    % unm1 = zeros(Ngx,Ngy);
    un   = zeros(Ngx,Ngy);
    unp1 = zeros(Ngx,Ngy);
    
    t = 0;
    % Update u0 and u1 to get started
    unm1 = u0(x,y);
%     % Approximation for u1
%     un(I,J) = unm1(I,J)+dt*u1(x(I,J),y(I,J))+(dt^2./(2*RHO(x(I,J),y(I,J))))...
%         .*(xDef(unm1,M11,I,J)+xyDef(unm1,M12,I,J)...
%         +yxDef(unm1,M21,I,J)+yDef(unm1,M22,I,J)+f(x(I,J),y(I,J),0));

    %New Approximation for u1
    % Approximation for u1
    un(I,J) = unm1(I,J)+dt*u1(x(I,J),y(I,J))+(dt^2./(2*rho(I,J)))...
        .*(xDef(unm1,m11,I,J)+xyDef(unm1,m12,I,J)...
        +yxDef(unm1,m21,I,J)+yDef(unm1,m22,I,J)+f(x(I,J),y(I,J),0));
    % Update boundary
    un(iax,I2b) = gax(y(iax,I2b),dt);
    un(ibx,I2b) = gbx(y(ibx,I2b),dt);
    un(I1b,iay) = gay(x(I1b,iay),dt);
    un(I1b,iby) = gby(x(I1b,iby),dt);
    % Check the approximation
    if(debug>1)
        uExact = ue(x,y,dt);
        err = uExact-un;
        errMax = max(max(abs(err)));
        fprintf('Debug: Check approx for U1: maxErr=%8.2e\n', errMax);
        return;
    end    
    


    %-----------------------Start time-stepping loop----------------------%
    for n = 1:Nt-1
        t = n*dt;
%         % Interior points
%         unp1(I,J) = 2*un(I,J)-unm1(I,J)+(dt^2./RHO(x(I,J),y(I,J))).*...
%             (xDef(un,M11,I,J)+xyDef(un,M12,I,J)+yxDef(un,M21,I,J)...
%             +yDef(un,M22,I,J)+f(x(I,J),y(I,J),t));
        
        % New interior points
        % Interior points
        unp1(I,J) = 2*un(I,J)-unm1(I,J)+(dt^2./rho(I,J)).*...
            (xDef(un,m11,I,J)+xyDef(un,m12,I,J)+yxDef(un,m21,I,J)...
            +yDef(un,m22,I,J)+f(x(I,J),y(I,J),t));
        % Update boundary
        unp1(iax,I2b) = gax(y(iax,I2b),t+dt);
        unp1(ibx,I2b) = gbx(y(ibx,I2b),t+dt);
        unp1(I1b,iay) = gay(x(I1b,iay),t+dt);
        unp1(I1b,iby) = gby(x(I1b,iby),t+dt);
        
        if(debug>0 && n==1)
            uExact = ue(x,y,t+dt);
            err = uExact-unp1;
            errMax = max(max(abs(err)));
            fprintf('Debug: after frist time-step: maxErr=%8.2e\n',errMax);
            return;
        end
        
       % Make a movie
       if (strcmp(coeffs,"const"))
           if (strcmp(ms, "trig"))
               step = 5;
           elseif (strcmp(ms, "pulse"))
               step = 2;
           elseif (strcmp(ms, "gauss"))
               step = 1;
           elseif (strcmp(ms, "gaussForce"))
               step = 1;
           else
               step = 2;
           end
       else
           if (strcmp(ms, "trig"))
               step = 5;
           else
               step = 1;
           end
       end
       if (makeMovie>0 && m==startInd && mod(n,step)==0)
           figure(1);
           surf(x,y,unp1);
           %surf(x,y,ue(x,y,t));
           colormap(getRainbow);
           xlabel("x"); ylabel("y"); zlabel("u");
           % shading interp;
           % view(0,90); % top view
           az = min(min(unp1));
           bz = max(max(unp1));
           title(sprintf("2D-Anisotropic Wave Equation t=%4.2e",t + dt));
           %axis([ax bx ay by az bz]);
           if (strcmp(coeffs,"const"))
               if (strcmp(ms,"trig"))
                   axis([ax bx ay by -3.5 3.5]);
               elseif (strcmp(ms, "pulse"))
                   axis([ax bx ay by -.5 5]);
               elseif (strcmp(ms, "trigPulse"))
                   axis([ax bx ay by -3 3]);
               elseif (strcmp(ms, "gauss"))
                   axis([ax bx ay by -1.5 1.5]);
               else
                   axis([ax bx ay by -.7 .7]);
               end
           else
               if (strcmp(ms,"trig"))
                   axis([ax bx ay by -6 5]);
               elseif (strcmp(ms, "pulse"))
                   axis([ax bx ay by -8 8]);
               elseif (strcmp(ms, "trigPulse"))
                   axis([ax bx ay by -5 5]);
               elseif (strcmp(ms, "gauss"))
                   axis([ax bx ay by -1 1]);
               else
                   axis([ax bx ay by -0.5 0.5]);
               end
           end
           if (strcmp(topView,"on"))
               %shading interp; 
               view(0,90); % top view
           end
           drawnow
           % Save movie
           if (makeMovie>1)
               % Get the current frame as an image
               frameImage = getframe(gcf);
               % Write the current frame to the video file
               writeVideo(videoWriterObj, frameImage);
           end
       end
       
        % Update unm1 & un for next step
        unm1 = un;
        un = unp1;
        
    end
    %-----------------------End time-stepping loop------------------------%
    
    %% Check error
    if (makeMovie==0 && ~startsWith(ms, "gauss") && (strcmp(coeffs,"const")||strcmp(ms,"poly")))
        % Compute error
        t = tFinal;
        uexact = ue(x,y,t);
        err = uexact-un;
        errMax(m) = max(max(abs(err)));
        maxNorm = max(max(un));
        NxList(m) = Nx; NtList(m) = Nt; dtList(m) = dt; relErr(m) = errMax(m)/maxNorm;
        fprintf('%s: t=%4.2f Nx=%3d Ny=%3d Nt=%3d dt=%8.2e maxNorm=%8.2e maxRelErr=%8.2e',...
            ms,tFinal,Nx,Ny,Nt,dt,maxNorm,errMax(m)/maxNorm);
        if( m>1 )
            fprintf(' ratio=%4.2f\n',errMax(m-1)/errMax(m));
        else
            fprintf('\n');
        end
        
        %% Plot results optional 
        if((m==3) && plotSol && (~strcmp(ms,'poly'))) % Plot for ms difers than poly
            figure(1);
            surf(x,y,un);
            colormap(getRainbow);
            colorbar;
            title(sprintf("2D Anisotropic Wave Equation solution, Nx=%d",Nx));
            xlabel("x"); ylabel("y"); zlabel("u");
            figure(2);
            surf(x,y,err);
            colormap(getRainbow);
            colorbar;
            title(sprintf("Error, N_x=%d",Nx));
            xlabel("x"); ylabel("y"); zlabel("err");
        end
        
    end
end
if (makeMovie>1)
    % Close the video file
    close(videoWriterObj);
end

%% Baseline: save or compare
ok = true;
if (makeMovie==0 && ~startsWith(ms,"gauss") && (strcmp(coeffs,"const")||strcmp(ms,"poly")))
    refDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'reference');
    refFile = fullfile(refDir, sprintf('matlab_%s_%s.mat', coeffs, ms));
    if saveRef
        if ~exist(refDir,'dir'), mkdir(refDir); end
        save(refFile, 'ms','coeffs','NxList','NtList','dtList','errMax','relErr','tFinal','cfl');
        fprintf('Saved %s\n', refFile);
        % Same values as text (BASELINE format of the C++ codes), for runBaselineSerial.sh
        txtFile = strrep(refFile, '.mat', '.txt');
        fid = fopen(txtFile, 'w');
        for i = startInd:numResolutions
            fprintf(fid, 'BASELINE: coeffs=%s ms=%s nx=%d numSteps=%d dt=%.16e maxErr=%.16e maxRelErr=%.16e\n', ...
                    coeffs, ms, NxList(i), NtList(i), dtList(i), errMax(i), relErr(i));
        end
        fclose(fid);
        fprintf('Saved %s\n', txtFile);
    elseif exist(refFile,'file')
        ref = load(refFile);
        d   = max(abs(errMax - ref.errMax));
        ok  = isequal(NtList, ref.NtList) && d < 1e-12;
        if ok, s = 'PASS'; else, s = 'FAIL'; end
        fprintf('Reference %s: max|errMax-ref| = %8.2e  %s\n', refFile, d, s);
    else
        ok = false;
        fprintf('Reference %s not found (run with saveRef = true first)\n', refFile);
    end
end