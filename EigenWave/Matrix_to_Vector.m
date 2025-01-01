function [V] = Matrix_to_Vector(U)

[n,m] = size(U);
V = zeros(n*m,1);
eqn = @(i,j) i+n*(j-1);
for i = 1:n
    for j = 1:m
        ie = eqn(i,j);
        V(ie) = U(i,j);
    end
end
% if (strcmp(geometry,'2D'))
%     for iy = 1:par.Ngy
%         for ix = 1:par.Ngx
%             ie = par.eqn(ix,iy);
%             V(ie) = U(ix,iy);
%         end
%     end
% elseif (strcmp(geometry,'Annulus'))
%     for ir = 1:par.Ngr
%         for ith = 1:par.Ngth
%             ie = par.eqn(ir,ith);
%             V(ie) = U(ir,ith);
%         end
%     end
% end
end