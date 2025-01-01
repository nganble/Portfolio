function int = getIntegral(U,dt,m)
int = 0;
for j = 1:m
    if((j==1)||(j==m))
        int = int+(1/2)*U(j);
    else
        int = int + U(j);
    end
end
int = int*dt;
end