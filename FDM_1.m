function call = FDM_1(s0,k,vol,r,t,nprice)

%explicit difference method
%nprice denote the number of price setps
%the number of time steps will be determined by conditional stability later

max_s = 2*max(s0, k);
dst = max_s / nprice;
ntime = floor(1.1*t*vol*vol*nprice*nprice);
dt = t / ntime;

st = zeros(nprice+1,1);
st(1) = 0;

for i = 2:nprice+1
    st(i) = st(i-1) + dst;
end


payoff = max(st-k,0);

grid = zeros(nprice+1,ntime+1);
grid(:,ntime+1) = payoff;

for i = 1:ntime
    k = ntime + 1 - i;
    for j = 2:nprice
        alp = dt * ((vol*vol*st(j)*st(j))/(2*dst*dst) + r*st(j)/(2*dst));
        beta = dt * (1/dt - (vol*vol*st(j)*st(j))/(dst*dst) - r);
        gamma = dt * ((vol*vol*st(j)*st(j))/(2*dst*dst) - r*st(j)/(2*dst));
        grid(j,k) = alp*grid(j+1,k+1) + beta*grid(j,k+1) + gamma*grid(j-1,k+1);
%       this can be alternatuvely done by following method            
%       delta = (grid(j+1,k+1) - grid(j-1,k+1))/(2*dst);
%       gamma = (grid(j+1,k+1) + grid(j-1,k+1) - 2*grid(j,k+1))/(dst*dst);
%       theta = -0.5*vol*vol*st(j)*st(j)*gamma - r*st(j)*delta + r*grid(j,k+1);
%       grid(j,k) = grid(j,k+1) - dt*theta;
    end
    grid(1,k) = 0;
    grid(nprice+1, k) = 2*grid(nprice, k) - grid(nprice-1, k);
end
call = grid(1+floor(s0/dst),1);

end