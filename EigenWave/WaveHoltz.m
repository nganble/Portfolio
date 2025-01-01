function [U] = WaveHoltz(V0,V1,par,fun)

% If option='testWaveSolver', return the solution to wave equation.
% Otherwise return the solution to WaveHoltz solver

global omega T option scheme c numCount geometry cpu
cpu0 = cputime;
numCount = numCount+1;
if (strcmp(geometry,'1D'))
    % Solver in 1D
    % Simplifying notation
    I   = par.I;
    dt  = par.dt;
    x   = par.x;
    iax = par.iax;
    ibx = par.ibx;
    Lh  = @(u,I) par.Lh(u,I);
    Ng  = par.Ngx;
    Nt  = par.Nt;
    
    % Allocate space for solutions at three levels
    unm1 = zeros(Ng,1);
    un   = zeros(Ng,1);
    unp1 = zeros(Ng,1);
    
    % Sol at t = 0
    unm1(I) = V0(I); 
    unm1(iax) = fun.ga(0);
    unm1(ibx) = fun.gb(0);
    
    % Apply time filter at time t = 0
    U = (1/2)*(cos(omega*(1-1)*dt)-1/4)*unm1;
    
    % Sol at t = dt
    if ((strcmp(scheme,'Implicit'))&&(strcmp(option,'EigenWave')))
        %un = cos(omega*dt)*unm1;
        un = par.dA\unm1;
    else
        un(I) = V0(I)+dt*V1(I)+(dt^2/2)*(c^2*Lh(V0,I)+fun.f(x(I),0));
        un(iax) = fun.ga(dt); % BC at x = a
        un(ibx) = fun.gb(dt); % BC at x = b
    end
    
 
    %------------------------Start time-stepping loop---------------------%
    for n=1:par.Nt-1
        if(strcmp(scheme,'Explicit'))
            t = n*dt; % new time
            unp1(I) = 2*un(I)-unm1(I)+(dt^2)*(c^2*Lh(un,I)+fun.f(x(I),t));
            % Update boundary point
            unp1(iax)=fun.ga(t+dt);
            unp1(ibx)=fun.gb(t+dt);
        elseif(strcmp(scheme,'Implicit'))
            rhs = zeros(Ng,1);
            rhs(I) = 2*un(I)-unm1(I)+((c*dt)^2/2)*Lh(unm1,I)...
                +(dt^2/2)*(fun.f(x(I),(n+1)*dt)+fun.f(x(I),(n-1)*dt));
            rhs(iax) = fun.ga((n+1)*dt); % BC at x = a
            rhs(ibx) = fun.gb((n+1)*dt); % BC at x = b
            unp1 = par.dA\rhs;
            %         % Check error after one step
            %           uexact = fun.ue(x,(n+1)*dt);
            %           errMax = max(abs(uexact-unp1));
            %           fprintf("maxErr = %5.2e\n",errMax);
            
        end
        % Applying the time-filter at t = n*dt
        U = U+(cos(omega*n*dt)-1/4)*un;
        % Update unm1 and un for next step
        unm1 = un;
        un = unp1;
        
    end
    %------------------------End time-stepping loop-----------------------%
    
    
elseif (strcmp(geometry,'2D'))
    % Solver in 2D
    % Simplify notations
    Ngx = par.Ngx; Ngy = par.Ngy; Nt = par.Nt;
    I  = par.I;  J  = par.J;
    Ib = par.Ib; Jb = par.Jb;
    iax = par.iax; ibx = par.ibx;
    iay = par.iay; iby = par.iby;
    x = par.x; y = par.y;
    dt = par.dt;
    
    % Allocate space for solutions at three levels
    unm1 = zeros(Ngx,Ngy);
    un   = zeros(Ngx,Ngy);
    unp1 = zeros(Ngx,Ngy);
    
    % Sol at t = 0
    unm1(I,J) = V0(I,J);
    unm1(iax,Jb) = fun.gax(y(iax,Jb),0);
    unm1(ibx,Jb) = fun.gbx(y(ibx,Jb),0);
    unm1(Ib,iay) = fun.gay(x(Ib,iay),0);
    unm1(Ib,iby) = fun.gby(x(Ib,iby),0);
    
    % Apply time filter at time t=0
    U = (1/2)*(cos(omega*(1-1)*dt)-1/4)*unm1;
    
    % Sol at t = dt
    if ((strcmp(scheme,'Implicit'))&&(strcmp(option,'EigenWave')))
        %un = cos(omega*dt)*unm1;
        w = par.dA\(Matrix_to_Vector(unm1));
        un = Vector_to_Matrix(w,par);
    else
        un(I,J) = unm1(par.I,par.J)+dt*V1(I,J)+...
            ((c*dt)^2/2)*par.Lh(unm1,I,J)+(dt^2)/2*fun.f(x(I,J),y(I,J),0);
        un(iax,Jb) = fun.gax(y(iax,Jb),dt);
        un(ibx,Jb) = fun.gbx(y(ibx,Jb),dt);
        un(Ib,iay) = fun.gay(x(Ib,iay),dt);
        un(Ib,iby) = fun.gby(x(Ib,iby),dt);
    end

    %     % Check approximation for first time-stepping
    %         uexact = fun.ue(x,y,dt);
    %         maxErr = max(max(abs(uexact-un)));
    %         fprintf('Check first step: errMax = %5.2e\n',maxErr);
    %
    
    %-------------------Start time-stepping loop----------------------%
    for n = 1:Nt-1
        t = n*dt; % New time
        if (strcmp(scheme,'Explicit'))
            unp1(I,J) = 2*un(I,J)-unm1(I,J)+(c*dt)^2*par.Lh(un,I,J)...
                +(dt^2)*fun.f(x(I,J),y(I,J),t);
            % Update boundary
            unp1(iax,Jb) = fun.gax(y(iax,Jb),t+dt);
            unp1(ibx,Jb) = fun.gbx(y(ibx,Jb),t+dt);
            unp1(Ib,iay) = fun.gay(x(Ib,iay),t+dt);
            unp1(Ib,iby) = fun.gby(x(Ib,iby),t+dt);
        elseif (strcmp(scheme,'Implicit'))
            rhs = zeros(Ngx,Ngy);
            rhs(I,J) = 2*un(I,J)-unm1(I,J) + 0.5*(dt^2)*par.Lh(unm1,I,J)...
                +(dt^2)/2*(fun.f(x(I,J),y(I,J),t-dt)+fun.f(x(I,J),y(I,J),t+dt));
            % Update boundary for rhs
            rhs(iax,Jb) = fun.gax(y(iax,Jb),t+dt);
            rhs(ibx,Jb) = fun.gbx(y(ibx,Jb),t+dt);
            rhs(Ib,iay) = fun.gay(x(Ib,iay),t+dt);
            rhs(Ib,iby) = fun.gby(x(Ib,iby),t+dt);
            
            rhs  = Matrix_to_Vector(rhs);
            unp1 = par.dA\rhs;
            unp1 = Vector_to_Matrix(unp1,par);
        end
        %             % Check error after ome step
        %             uexact = fun.ue(x,y,t+dt);
        %             maxErr = max(max(abs(uexact-unp1)));
        %             fprintf("errMax = %5.2e\n",maxErr);
        
        % Applying the time-filter at t=n*dt
        U = U+(cos(omega*n*dt)-1/4)*un;
        % Update umn1 & un for next step
        unm1 = un;
        un   = unp1;
        
    end
    %----------------------End time-stepping loop---------------------%
    
    
elseif (strcmp(geometry,'Annulus'))
    % Solver on the annulus
    % Simplify notations
    I    = par.I;    Ib   = par.Ib;
    JbL  = par.JbL;  Jb   = par.Jb; J = par.J;
    iar  = par.iar;  ibr  = par.ibr;
    iath = par.iath; ibth = par.ibth;
    iaGhost = par.iaGhost;
    Ngr = par.Ngr; Ngth = par.Ngth; Nt = par.Nt;
    dt = par.dt;
    th = par.th;
    x  = par.x;
    y  = par.y;
    Lh = @(u,I,J) par.Lh(u,I,J);
    
    % Allocate space for solutions of three levels
    unm1 = zeros(Ngr,Ngth);
    un   = zeros(Ngr,Ngth);
    unp1 = zeros(Ngr,Ngth);
    
    % Sol at t = 0
    unm1(I  ,JbL   ) = V0(I,JbL);
    unm1(iar,Jb    ) = fun.ga(th(iar,Jb),0);
    unm1(ibr,Jb    ) = fun.gb(th(ibr,Jb),0);
    unm1(Ib,iaGhost) = unm1(Ib,ibth-1);
    unm1(Ib,ibth   ) = unm1(Ib,iath);
    
    % Apply time filter at time t = 0
    U = (1/2)*(cos(omega*(1-1)*dt)-1/4)*unm1;
    
    % Sol at t = dt
    if ((strcmp(scheme,'Implicit'))&&(strcmp(option,'EigenWave')))
        %un(I,JbL) = cos(omega*dt)*unm1(I,JbL);
        rhs = zeros(Ngr,Ngth);
        rhs(I,J) = unm1(I,J);
        rhs(iar,Jb    ) = fun.ga(th(iar,Jb),dt);
        rhs(ibr,Jb    ) = fun.gb(th(ibr,Jb),dt);
        rhs(Ib,iaGhost) = zeros(Ngr,1);
        rhs(Ib,ibth   ) = zeros(Ngr,1);
        rhs = Matrix_to_Vector(rhs);
        w = par.dA\rhs;
        un = Vector_to_Matrix(w,par);
    else
        un(I,JbL) = unm1(I,JbL)+dt*V1(I,JbL)+((c*dt)^2/2)*Lh(unm1,I,JbL)...
            +dt^2/2*fun.f(x(I,JbL),y(I,JbL),0);
        un(iar,Jb    ) = fun.ga(th(iar,Jb),dt);
        un(ibr,Jb    ) = fun.gb(th(ibr,Jb),dt);
        un(Ib,iaGhost) = un(Ib,ibth-1);
        un(Ib,ibth   ) = un(Ib,iath);
    end
    
    %     % Check approximation at t = dt
    %     if (strcmp(ms,'eigen'))
    %         uexact = fun.ue(r,th,dt);
    %     else
    %         uexact = fun.ue(x,y,dt);
    %     end
    %     maxErr = max(max(abs(uexact-un)));
    %     fprintf('Check first step: errMax = %5.2e\n',maxErr);
    %     pause;
    
    % ------------------ Start time-stepping loop ------------------- %
    for n = 1:Nt-1
        t = n*dt;
        if (strcmp(scheme,'Explicit'))
            unp1(I,JbL) = 2*un(I,JbL)-unm1(I,JbL)+((c*dt)^2)*Lh(un,I,JbL)...
                +dt^2*fun.f(x(I,JbL),y(I,JbL),t);
            unp1(iar,Jb   ) = fun.ga(th(iar,Jb),dt+t);
            unp1(ibr,Jb   ) = fun.gb(th(ibr,Jb),dt+t);
            unp1(Ib,iath-1) = unp1(Ib,ibth-1);
            unp1(Ib,ibth  ) = unp1(Ib,iath);
        elseif (strcmp(scheme,'Implicit'))
            rhs = zeros(Ngr,Ngth);
            rhs(I,JbL) = 2*un(I,JbL)-unm1(I,JbL)+((c*dt)^2/2)*Lh(unm1,I,JbL)...
                +(dt^2)/2*(fun.f(x(I,JbL),y(I,JbL),t-dt)+fun.f(x(I,JbL),y(I,JbL),t+dt));
            % Update boundary for rhs
            rhs(iar,Jb    ) = fun.ga(th(iar,Jb),t+dt);
            rhs(ibr,Jb    ) = fun.gb(th(ibr,Jb),t+dt);
            rhs(Ib,iaGhost) = zeros(Ngr,1);
            rhs(Ib,ibth   ) = zeros(Ngr,1);
            
            rhs  = Matrix_to_Vector(rhs);
            unp1 = par.dA\rhs;
            unp1 = Vector_to_Matrix(unp1,par);
        end
        
        % Applying time filter at t=n*dt
        U = U+(cos(omega*n*dt)-1/4)*un;
        
        % Update for next step
        unm1 = un;
        un   = unp1;
    end
    %------------------ End time-stepping loop ---------------------- %
    
end

    % Apply time-filter at final time step
    U = U+(1/2)*(cos(omega*Nt*dt)-1/4)*unp1;
    U = (2/T)*dt*U;
    
    
if(strcmp(option,'testWaveSolver'))
    U = unp1;
end
cpu = cputime-cpu0;

end