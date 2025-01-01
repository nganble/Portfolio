function [a] = dis_norm(U,par)
a = dis_product(U,U,par);
a = sqrt(a);
end