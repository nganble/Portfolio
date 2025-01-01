function fun = getMS(dm,ms)
global option c geometry
if(strcmp(option,'testWaveSolver'))
    %% MS in 1D
    if(strcmp(geometry,'1D'))
        if(strcmp(ms,'poly'))
            % Polynomial ms
            b0 = .5; b1 = .7; b2 = .9; % time coeffivelcients
            c0 = 1;  c1 = .8; c2 = .6;  % space coefficients
            fun.ue = @(x,t) (b0+b1*t+b2*t^2)*(c0+c1*x+c2*x.^2);
            uet    = @(x,t) (b1+2*b2*t)*(c0+c1*x+c2*x.^2);
            uett   = @(x,t) 0*t+(2*b2)*(c0+c1*x+c2*x.^2);
            uexx   = @(x,t) (b0+b1*t+b2*t^2).*2*c2+0*x;
        elseif(strcmp(ms,'pulse'))
            beta = 20; x0 = .25;
            fun.ue = @(x,t) exp(-(beta^2)*(x-x0-c*t).^2);
            uet    = @(x,t) 2*c*beta^2*(x-x0-c*t).*fun.ue(x,t);
            uett   = @(x,t) (4*c^2*beta^4*(x-x0-c*t).^2-2*c^2*beta^2).*fun.ue(x,t);
            uexx   = @(x,t) (4*beta^4*(x-x0-c*t).^2-2*beta^2).*fun.ue(x,t);
        elseif(strcmp(ms,'eigen'))
            % Test Eigenfunction
            n = 1;
            lambda = n*pi;
            fun.ue = @(x,t) cos(lambda*c*t)*sin(lambda*x);
            uet    = @(x,t) -lambda*c*sin(lambda*c*t)*sin(lambda*x);
            uett   = @(x,t) -lambda^2*c^2*cos(lambda*c*t)*sin(lambda*x);
            uexx   = @(x,t) -lambda^2*cos(lambda*c*t)*sin(lambda*x);
        else
            %         fun.ue = @(x,t) sin(x)*(t^3);
            %         uet = @(x,t) 3*sin(x)*t^2;
            %         uett = @(x,t) 6*sin(x)*t;
            %         uexx = @(x,t) -sin(x)*t^3;
            fun.ue = @(x,t) sin(x+c*t)+cos(x-c*t);
            uet    = @(x,t) cos(x+c*t)-sin(x-c^t);
            uett   = @(x,t) 0*x+0*t;
            uexx   = @(x,t) 0*x+0*t;
            
        end
        % Forcing term
        fun.f = @(x,t) uett(x,t)-c^2*uexx(x,t);
        % Boundary condition
        fun.ga = @(t) fun.ue(dm.LB,t);
        fun.gb = @(t) fun.ue(dm.RB,t);
        % Initial condition
        fun.u0 = @(x) fun.ue(x,0);
        fun.u1 = @(x) uet(x,0);
    elseif(strcmp(geometry,'2D'))
        %% MS in 2D
        if(strcmp(ms,'poly'))
            b0  = .5;  b1 = .7;  b2  = .9; % time coefficients
            c00 =  1; c10 = .8;  c20 = .6; % space coefficients in x
            c01 = .3; c11 = .25; c02 = .2; % space coefficients in y
            % Polynomial in t
            polyT    = @(x,y,t) b0+b1*t+b2*t^2+0*x+0*y;
            dtPolyT  = @(x,y,t) b1+2*b2*t+0*x+0*y;
            dttPolyT = @(x,y,t) 2*b2+0*x+0*y;
            % Polynomial in x&y
            polyXY    = @(x,y,t) c00+c10*x+c20*x.^2+c11*x.*y+c01*y+c02*y.^2;
            dxxPolyXY = @(x,y,t) 2*c20+0*x+0*y;
            dyyPolyXY = @(x,y,t) 2*c02+0*x+0*y;
            % Uexact and its derivaties
            fun.ue = @(x,y,t) polyT(x,y,t).*polyXY(x,y,t);
            uet    = @(x,y,t) dtPolyT(x,y,t).*polyXY(x,y,t);
            uett   = @(x,y,t) dttPolyT(x,y,t).*polyXY(x,y,t);
            uexx   = @(x,y,t) polyT(x,y,t).*dxxPolyXY(x,y,t);
            ueyy   = @(x,y,t) polyT(x,y,t).*dyyPolyXY(x,y,t);
        elseif(strcmp(ms,"trig"))
            % Trigonometric manufactured solution
            kx = 3*pi;
            ky = 2*pi;
            gamma = c*sqrt(kx^2+ky^2);
            % Function in t
            funT    = @(t) sin(gamma*t);
            dtFunT  = @(t) gamma*cos(gamma*t);
            dttFunT = @(t) -gamma^2*sin(gamma*t);
            % Function in x
            funX    = @(x) sin(kx*x);
            dxxFunX = @(x) -kx^2*sin(kx*x);
            % Function in y
            funY    = @(y) sin(ky*y);
            dyyFunY = @(y) -ky^2*sin(ky*y);
            % Exact solution
            fun.ue = @(x,y,t) funT(t)*funX(x).*funY(y);
            uet    = @(x,y,t) dtFunT(t)*funX(x).*funY(y);
            uett   = @(x,y,t) dttFunT(t)*funX(x).*funY(y);
            uexx   = @(x,y,t) funT(t)*dxxFunX(x).*funY(y);
            ueyy   = @(x,y,t) funT(t)*funX(x).*dyyFunY(y);
        else % Eigenfunction
            lambda_n = 2;
            lambda_m = 2;
            gamma = c*sqrt((lambda_n*pi)^2+(lambda_m*pi)^2);
            fun.ue = @(x,y,t) cos(gamma*t)*sin(lambda_n*pi*x).*sin(lambda_m*pi*y);
            uet    = @(x,y,t) -gamma*sin(gamma*t)*sin(lambda_n*pi*x).*sin(lambda_m*pi*y);
            uett   = @(x,y,t) -gamma^2*fun.ue(x,y,t);
            uexx   = @(x,y,t) -(lambda_n*pi)^2*fun.ue(x,y,t);
            ueyy   = @(x,y,t) -(lambda_m*pi)^2*fun.ue(x,y,t);
        end
        
        % Initial conditions
        fun.u0 = @(x,y) fun.ue( x,y,0);
        fun.u1 = @(x,y) uet(x,y,0);
        % Boundary conditions
        fun.gax = @(y,t) fun.ue(dm.LB,y,t);
        fun.gbx = @(y,t) fun.ue(dm.RB,y,t);
        fun.gay = @(x,t) fun.ue(x,dm.BB,t);
        fun.gby = @(x,t) fun.ue(x,dm.TB,t);
        % Forcing term
        fun.f = @(x,y,t) uett(x,y,t)-c^2*(uexx(x,y,t)+ueyy(x,y,t));
    elseif(strcmp(geometry,'Annulus'))
        %% MS on the annulus
        if(strcmp(ms,'poly'))
            % Polynomial manufactured solution
            b0  = .5; b1  = .7;  b2  = .9; % time coefficients
            c00 = 1;  c10 = .8;  c20 = .6; % space coefficients
            c01 = .3; c11 = .25; c02 = .2;
            % Polynomial in t
            polyT    = @(x,y,t) b0+b1*t+b2*t^2+0*x+0*y;
            dTpolyT  = @(x,y,t) b1+2*b2*t+0*x+0*y;
            dTTpolyT = @(x,y,t) 2*b2+0*x+0*y;
            % Polynomial in x&y
            polyXY    = @(x,y,t) c00+c10*x+c20*x.^2+c11*x.*y+c01*y+c02*y.^2;
            dXXpolyXY = @(x,y,t) 2*c20+0*x+0*y;
            dYYpolyXY = @(x,y,t) 2*c02+0*x+0*y;
            % Uexact and its derivaties
            fun.ue = @(x,y,t) polyT(x,y,t).*polyXY(x,y,t);
            uet    = @(x,y,t) dTpolyT(x,y,t).*polyXY(x,y,t);
            uett   = @(x,y,t) dTTpolyT(x,y,t).*polyXY(x,y,t);
            uexx   = @(x,y,t) polyT(x,y,t).*dXXpolyXY(x,y,t);
            ueyy   = @(x,y,t) polyT(x,y,t).*dYYpolyXY(x,y,t);
            
        elseif(strcmp(ms,"trig"))
            %% Trigonometric manufactured solution
            kx = pi;
            ky = pi;
            gamma = c*sqrt(kx^2+ky^2);
            % Function of t
            funT    = @(t) sin(gamma*t);
            dTfunT  = @(t) gamma*cos(gamma*t);
            dTTfunT = @(t) -gamma^2*sin(gamma*t);
            % Function of x
            funX    = @(x) sin(kx*x);
            dXXfunX = @(x) -kx^2*sin(kx*x);
            % Function of y
            funY    = @(y) sin(ky*y);
            dYYfunY = @(y) -ky^2*sin(ky*y);
            % Uexact and its derivatives
            fun.ue = @(x,y,t) funT(t)*funX(x).*funY(y);
            uet    = @(x,y,t) dTfunT(t)*funX(x).*funY(y);
            uett   = @(x,y,t) dTTfunT(t)*funX(x).*funY(y);
            uexx   = @(x,y,t) funT(t)*dXXfunX(x).*funY(y);
            ueyy   = @(x,y,t) funT(t)*funX(x).*dYYfunY(y);
        elseif (strcmp(ms,"eigen"))
            %% Eigenfunction
            m = 1;
            %lambda = 2.2237527968057966;
            lambda = 2.2217160733567396;
            phi = @(u,v) (besselj(m,lambda*u)-besselj(m,lambda*dm.LB)...
                /bessely(m,lambda*dm.LB)*bessely(m,lambda*u)).*cos(m*v);
            fun.ue = @(r,th,t) phi(r,th).*cos(lambda*t);
            uet    = @(r,th,t) -lambda*phi(r,th).*sin(lambda*t);
        end
        polarUe = @(r,theta,t) fun.ue(r.*cos(theta),r.*sin(theta),t);
        % Initial conditions
        fun.u0 = @(x,y) fun.ue(x,y,0);
        fun.u1 = @(x,y) uet(x,y,0);
        % Boundary conditions and forcing
        if (strcmp(ms,"eigen"))
            fun.ga = @(theta,t) fun.ue(dm.LB,theta,t);
            fun.gb = @(theta,t) fun.ue(dm.RB,theta,t);
            fun.f  = @(x,y,t) 0*x+0*y+0*t;
        else
            fun.ga = @(theta,t) polarUe(dm.LB,theta,t);
            fun.gb = @(theta,t) polarUe(dm.RB,theta,t);
            fun.f  = @(x,y,t) uett(x,y,t)-c^2*(uexx(x,y,t)+ueyy(x,y,t));
        end
    end
else
    %% Functions for WaveHoltz solver
    if (strcmp(geometry,'1D'))
        fun.ga = @(t) 0*t;
        fun.gb = @(t) 0*t;
        fun.f  = @(x,t) 0.*x+0*t;
    elseif (strcmp(geometry,'2D'))
        fun.gax = @(y,t) 0*y+0*t;
        fun.gbx = @(y,t) 0*y+0*t;
        fun.gay = @(x,t) 0*x+0*t;
        fun.gby = @(x,t) 0*x+0*t;
        fun.f   = @(x,y,t) 0*x+0*y+0*t;
    elseif (strcmp(geometry,'Annulus'))
        fun.ga = @(theta,t) theta*0+t*0;
        fun.gb = @(theta,t) theta*0+t*0;
        fun.f  = @(x,y,t) 0*x+0*y+0*t;
    end
    
    
end
