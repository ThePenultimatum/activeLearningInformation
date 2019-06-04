clear;
global T b x0 sigma mu epsilon resolution trajectory L dt phix K rowres colres;

T = 1000;
resolution = 1;
dt = 1/resolution;
b = 0;
bmax = 2;
numBsamples = 100;


%x0 = [0; 1];
x0 = [1; 2];

sigma = [1 0; 0 1];

mu = 1;%[1; 1];%0;
L = 2;%2;

%epsilon = 0.1;

K = 10;

rowres = 10;
colres = 10;

%trajectory = [0; 1];
%currTraj = x0;
%for i=1:T*resolution
%    next = currTraj + dt * xdot(currTraj);
%    trajectory(:,i) = next;
%    currTraj = next;
%end
% end
% [t, x] = ode45(@(t,x) [0 1; -1 -b]*x, linspace(0,100, 101), x0);
% trajectory = x;

%%%%

%i = 0;
%%%% 2 norm is sqrt of (parts squared)

gammaKs = []; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% is this a matrix instead of a vector?????????????
for ki=0:K
    for kj=0:K
        %newone = getGammaK(k);
        %gammaKs = [gammaKs; newone];
        newone = getGammaKiKj(ki, kj);
        gammaKs(ki+1,kj+1) = newone;
    end
end

gammaIJs = [];
hs = [];
fkxs = [];
cks = [];
phiks = [];

phix = phi();%(trajectory);

%for t=1:((T/dt) + 1)
%    t
%    xs = trajectory(:,t);

bsToTry = linspace(0,bmax,numBsamples);
bsTried = [];
epsilonsRes = [];
epsilonProgress = [];

for i=1:(length(bsToTry))
    
    b = bsToTry(i);
    bsTried(i) = b;
    
    
    [t, x] = ode45(@(t,x) [0 1; -1 -b]*x, linspace(0,T,T+1), x0);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    trajectory = x;
    
    sumSoFar = 0;
    
    for xk=0:K
        for yk=0:K
            ks = [xk; yk];
            %
            gammaIJ = gammaKs(xk+1, yk+1);
            gammaIJs(xk+1, yk+1) = gammaIJ;
            %
            h = getHK(ks);
            hs(xk+1, yk+1) = h;
            %
            fkx = getFkx(trajectory, ks, h);
            fkxs = [fkxs ; fkx];
            %
            ck = getCks(fkx);
            cks(xk+1, yk+1) = ck;
            %
            phik = getPhik(ks, h);
            phiks(xk+1, yk+1) = phik;
            %
            %phik = getPhik(fkx);
            sumSoFar = sumSoFar + gammaIJ * (norm(ck - phik))^2;
            epsilonProgress(xk+1, yk+1) = sumSoFar;
        end
    end
    
    epsilon_t = sumSoFar;
    epsilonsRes(i) = epsilon_t;
end
%end


plot(bsTried, epsilonsRes);
title("Epsilon vs b");
xlim([0 bmax]);
xlabel("b");
ylabel("epsilon");

%%%% functions

function pkx = getPhik(ks, hk)%(fkx)
    global phix rowres colres L;
    %phix;
    %fkx;
    %pkx = sum(phix) * transpose(fkx);
    
    sumsofar = 0;
    drow = L / rowres;
    dcol = L / colres;
    for a=0:(rowres+1)
        for b=0:(colres+1)
            drowab = drow * a;
            dcolab = dcol * b;
            elem = (1/hk)*cos(ks(1)*pi*drowab/L) * cos(ks(2)*pi*dcolab/L);
            sumsofar = sumsofar + elem;
        end
    end
    
    pkx = phix(ks(1)+1, ks(2)+1)*sumsofar;%/(colres*rowres);
end

function ck = getCks(fkx)
    global T dt;
    fkx;
    %fun = @(t) fkx(t);
    %ck = integral(fun, 1, T/dt+1); %0;%(1/T) *
    %ck = (1/T) * [sum(fkx(1,:)) ; sum(fkx(2,:))];
    ck = (1/T) * sum(fkx);%(1,:));
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
    x2part = cos(k(2) * pi * x(2,:) / L);
    tot = resi * (x1part .* x2part);
    fkx = tot;
end

function hk = getHK(ks)
    global L 
    fun = @(x1,x2) (cos(ks(1) * pi * x1 / L).^2).*((cos(ks(2) * pi * x2 / L)).^2);
    hk = sqrt(integral2(fun, 0, L, 0, L));
end

function g = getGammaKiKj(i, j) 
    %%%%%%%%%%% this assumes that k is a column vector input or a scalar
    usek = [i; j];
    n = 2;
    g = (1+norm(usek)^2)^(-1*((n+1)/2));
end

function phi_x = phi()
  global sigma mu L K;
  newk = K + 1;
  xval = [linspace(0, L, newk); linspace(0, L, newk)];
  phi_x = (1/(sqrt(det(2*pi*sigma))))*exp(-0.5*transpose(xval-mu)*inv(sigma)*(xval-mu));
end