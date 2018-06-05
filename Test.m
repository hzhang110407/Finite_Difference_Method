%set parameters
s0 = 100; %today's price
k = 100; %strike
t = 1; %time to expiry
vol = 0.2; %volatility
r = 0.05; % risk-free rate

%get the closed form solution
call = bs(s0, k, vol, t, 1);
%run each version 
call_1 = FMD_1(s0, k, vlo. r, t, 100);
call_2 = FDM_2(s0, k, vol, r, t, 100, 100);
call_3 = FDM_3(s0, k, vol, r, t, 100, 100, 0.5);
call_4 = FDM_4(s0, k, vol, r, t, 100, 100, 0.5);


function [call, put] = bs(stock, strike, rate, time, volatility)

    d1 = (1/(volatility*power(time, 0.5))) * ( -log(strike/stock) + (rate + power(volatility,2)/2) * time);
    d2 = d1 - volatility*power(time, 0.5);

    call = normcdf(d1)*stock - normcdf(d2)*exp(-rate*time).*strike;
    put = strike*exp(-rate*time) - stock + call;

end

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

function call = FDM_2(s0,k,vol,r,t,nprice,ntime)

    %implicit difference method

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
        target = zeros(nprice+1,1);
        target(2:nprice) = grid(2:nprice,k+1);

        v1 = dt/(dst*dst);
        v2 = dt/dst;

        m = zeros(nprice+1);
        m(nprice+1,nprice-1:nprice+1) = [1, -2, 1];%to be modified, when s is large, V is almost linear respect to s
        m(1, 1) = 1; %to be modified, V(0,t) = 0 for any T >= t >= 0
        for j = (2:nprice)
            m(j,j-1) = -0.5*vol*vol*st(j)*st(j)*v1 + 0.5*r*st(j)*v2;
            m(j,j) = 1 + dt*r + vol*vol*st(j)*st(j)*v1;
            m(j,j+1) = -0.5*vol*vol*st(j)*st(j)*v1 - 0.5*r*st(j)*v2;
        end
        %disp(m);
        grid(:,k) = linsolve(m, target);

    end
    call = grid(1+floor(s0/dst),1);    
end

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

function call = FDM_4(s0,k,vol,r,t,nprice,ntime,w)

    call_1 = FDM_3(s0,k,vol,r,t,nprice,ntime,w);
    call_2 = FDM_3(s0,k,vol,r,t,2*nprice,4*ntime,w);
    call = (4*call_2 - call_1)/3;

end
