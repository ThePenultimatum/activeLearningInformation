global T b x0 sigma mu epsilon;

T = 100;
b = 0;

x0 = [0; 1];

sigma = Diag(2,2);

mu = 0;

epsilon = 0.1;

%%%%

i = 0;
%%%% 2 norm is sqrt of (parts squared)
while sqrt((dirDerivJ(eta_i, zeta_i))^2) < epsilon
    i = i + 1;
end

%%%% functions

function dj = dirDerivJ(eta_i, zeta_i)
  dj = 0;
end

function xd = xdot(x)
  global b;
  xd = [0 1; -1 -b] * x;
end

function phi_x = phi(x)
  muhere = [0; 0]; % same shape as x
  phi_x = (1/(sqrt(det(2*pi*sigma))))*exp(-0.5*transpose(x-muhere)*inv(sigma)*(x-muhere));
end