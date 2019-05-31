global T b x0 sigma mu epsilon resolution trajectory L dt;

T = 100;
resolution = 1;
dt = 1/resolution;
b = 0;

x0 = [0; 1];

sigma = [1 0; 0 1];

mu = 0;
L = 2;

epsilon = 0.1;

K = 10;

trajectory = [0; 1];
currTraj = x0;
for i=1:T*resolution
    next = currTraj + dt * xdot(currTraj);
    trajectory(:,i) = next;
    currTraj = next;
end

%%%%

%i = 0;
%%%% 2 norm is sqrt of (parts squared)

gammaKs = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% is this a matrix instead of a vector?????????????
for ki=1:K
    for kj=1:K
        %newone = getGammaK(k);
        %gammaKs = [gammaKs; newone];
        newone = getGammaKiKj(ki, kj);
        gammaKs(ki,kj) = newone;
    end
end

%for t=1:((T/dt) + 1)
%    t
%    xs = trajectory(:,t);
    for xk=1:K
        for yk=1:K
            ks = [xk; yk];
            gammaIJ = gammaKs(xk, yk);
            h = getHK(ks);
            fkx = getFkx(trajectory, ks, h)
            ck = getCks(fkx);
            phik = [];
        end
    end
%end

%%%% functions

function ck = getCks(fkx)
    global T dt;
    fkx
    %fun = @(t) fkx(t);
    %ck = integral(fun, 1, T/dt+1); %0;%(1/T) *
    ck = (1/T) * sum(fkx(1,:));
end

function fkx = getFkx(x, k, h)
    global L
    1;
    %h = getHK(k);
    %L = 1;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%????????????????????????????/
    res = [];
    normalizer = 1/h;
    resi = normalizer * 1;
    
    x1part = cos(k(1) * pi * x(1,:) / L);
    x2part = cos(k(1) * pi * x(2,:) / L);
    tot = resi * (x1part .* x2part);
    %for i=1:2
    %    ki = k(i);
    %    resi = resi * cos(k(i) * pi * x(:,i) / L);
    fkx = tot;
    %end
    
    %fkx = 0;
end

function hk = getHK(ks)
    global L 
    fun = @(x1,x2) (cos(ks(1) * pi * x1 / L).^2).*((cos(ks(2) * pi * x2 / L)).^2);
    hk = sqrt(integral2(fun, 0, 2, 0, 2));
end

function g = getGammaKiKj(i, j) 
    %%%%%%%%%%% this assumes that k is a column vector input or a scalar
    usek = [i; j];
    n = 2;
    g = (1+sqrt(transpose(usek)*usek)^2)^(-1*((n+1)/2));
end

function xd = xdot(x)
  global b;
  xd = [0 1; -1 -b] * x;
end

function phi_x = phi(x)
  global sigma mu;
  muhere = [0; 0]; % same shape as x
  phi_x = (1/(sqrt(det(2*pi*sigma))))*exp(-0.5*transpose(x-mu)*inv(sigma)*(x-mu));
end