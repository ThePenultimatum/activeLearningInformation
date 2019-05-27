global T b x0 sigma mu epsilon;

T = 100;
b = 0;

x0 = [0; 1];

sigma = Diag(2,2);

mu = 0;

epsilon = 0.1;

%%%%

i = 0;

while ((dirDeriv(J, eta_i, zeta_i))^2) < epsilon
    i++;
end

%%%% functions

function xd = xdot(x)
  global b;
  xd = [0 1; -1 -b] * x;
end

function phi_x = phi(x)
  muhere = [0; 0]; % same shape as x
  phi_x = (1/(sqrt(det(2*pi*sigma))))*exp(-0.5*transpose(x-muhere)*inv(sigma)*(x-muhere));
end