function call = FDM_3(s0,k,vol,r,t,nprice,ntime,w)

dt = t/ntime; %this method is unconditionally stable, no need to adjust time step accordingly
max_s = 2*max(s0, k);
dst = max_s / nprice;

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
    
    v1 = dt/(dst*dst);
    v2 = dt/dst;
    m = zeros(nprice-1, nprice+1);
    
    for j = 1:nprice-1
        A = vol*vol*st(j+1)*st(j+1);
        B = r*st(j+1);
        m(j,j) = 0.5*A*v1*(w-1) - 0.5*B*v2*(w-1);
        m(j,j+1) = -1 - A*v1*(w-1) - dt*r*(w-1);
        m(j,j+2) = 0.5*A*v1*(w-1) + 0.5*B*v2*(w-1);        
    end
    
    target = zeros(nprice+1,1);
    target(2:nprice) = m*grid(:,k+1);
    
    m = zeros(nprice+1);
    m(nprice+1,nprice-1:nprice+1) = [1, -2, 1];%to be modified, when s is large, V is almost linear respect to s
    m(1, 1) = 1; %to be modified, V(0,t) = 0 for any T >= t >= 0
    for j = (2:nprice)
        m(j,j-1) = 0.5*vol*vol*st(j)*st(j)*v1*w - 0.5*r*st(j)*v2*w;
        m(j,j) = -1 - dt*r*w - vol*vol*st(j)*st(j)*v1*w;
        m(j,j+1) = 0.5*vol*vol*st(j)*st(j)*v1*w + 0.5*r*st(j)*v2*w;
    end
    %disp(m);
    grid(:,k) = linsolve(m, target);
    
end

call = grid(1+floor(s0/dst),1);
    
end
