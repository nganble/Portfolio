function dm = getDomain(geometry)
% GETDOMAIN: domain bounds for each geometry
%   dm = getDomain(geometry)
%   geometry : '1D' | '2D' | 'Annulus'
%   dm       : LB, RB (x or r), and BB, TB (y or theta) in 2D/Annulus

if (strcmp(geometry,'1D'))
    dm.LB = 0;
    dm.RB = 1;
elseif(strcmp(geometry,'2D'))
    dm.LB = 0;
    dm.RB = 1;
    dm.BB = 0;
    dm.TB = 1;
elseif (strcmp(geometry,'Annulus'))
    dm.LB = 0.5;
    dm.RB = 2;
    dm.BB = 0;
    dm.TB = 2*pi;
end
end