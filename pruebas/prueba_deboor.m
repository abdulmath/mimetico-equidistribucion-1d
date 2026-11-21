% prueba_deboor.m
% Tres comprobaciones sobre el paso de de Boor y sobre el lazo exterior:
%   1. que un solo paso equidistribuya la densidad hasta el redondeo
%   2. cómo depende la convergencia del lazo del factor de relajación
%   3. la tabla de errores en la malla final para varios tamaños

function prueba_deboor()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

a = 0;  b = 1;  alpha = 250;  x0 = 0.8;
kf = @(z) 1/alpha + alpha*(z-x0).^2;
ff = @(z) -2*(1 + alpha*(z-x0).*(atan(alpha*x0) + atan(alpha*(z-x0))));
ue = @(z) (1-z).*(atan(alpha*(z-x0)) + atan(alpha*x0));
gu = @(z) -(atan(alpha*(z-x0)) + atan(alpha*x0)) ...
          + (1-z).*(alpha./(1+alpha^2*(z-x0).^2));
ga = alpha/(1+alpha^2*x0^2);
gb = -atan(alpha*(1-x0)) - atan(alpha*x0);

% --- 1. exactitud de un paso ---------------------------------------------
printf('\n=== Un paso de de Boor sobre una densidad conocida ===\n');
printf('  Se mide la dispersión relativa de la masa por celda, int rho sobre\n');
printf('  cada celda, que debe anularse hasta el redondeo.\n');
printf('  %-9s %-14s %-14s\n', 'Celdas', 'antes', 'después');
for N = [50 100 200 400]
  x  = linspace(a, b, N+1);
  xc = (x(1:N) + x(2:N+1))/2;
  r  = sqrt(1 + 0.5*gu(xc).^2);
  d0 = disper(masa(x, r, x));
  xn = deBoor(x, r);
  d1 = disper(masa(x, r, xn));
  printf('  %-9d %-14.3e %-14.3e\n', N, d0, d1);
end

% --- 2. efecto de la relajación ------------------------------------------
printf('\n=== Convergencia del lazo exterior según omega, 100 celdas ===\n');
printf('  %-8s %-9s %-8s %-13s %-13s\n', 'omega', 'converge', 'iter', 'mov. final', 'Linf');
x0m = linspace(a, b, 101);
for om = [1.0 0.7 0.5 0.3 0.2 0.1]
  o = struct('omega', om, 'maxIt', 300, 'tolX', 1e-10);
  [xa, Ua, h] = adaptaDeBoor(x0m, kf, ff, ga, gb, o);
  [eI, ~] = erroresMalla(xa, Ua, ue);
  if h.convergio, s = 'sí'; else, s = 'no'; end
  printf('  %-8.2f %-9s %-8d %-13.3e %-13.4e\n', om, s, h.it, h.mov(end), eI);
end

% --- 3. tabla final -------------------------------------------------------
printf('\n=== Errores en la malla final, omega = 0.3 ===\n');
printf('  %-8s %-7s %-7s %-13s %-13s %-13s %-13s %-8s\n', 'Celdas', 'nodos', ...
       'iter', 'Linf adapt', 'L2 adapt', 'Linf unif', 'L2 unif', 'ganancia');
o = struct('omega', 0.3, 'maxIt', 300, 'tolX', 1e-10);
for N = [40 60 80 100 120 160 200 400]
  xu = linspace(a, b, N+1);
  [xa, Ua, h] = adaptaDeBoor(xu, kf, ff, ga, gb, o);
  [eIa, e2a] = erroresMalla(xa, Ua, ue);
  Uu = resuelveRobin1D(xu, kf, ff, ga, gb);
  [eIu, e2u] = erroresMalla(xu, Uu, ue);
  printf('  %-8d %-7d %-7d %-13.4e %-13.4e %-13.4e %-13.4e %-8.1f\n', ...
         N, N+1, h.it, eIa, e2a, eIu, e2u, eIu/eIa);
end
printf('\n');
end

% --------------------------------------------------------------------------
function d = disper(v)
  d = (max(v) - min(v))/mean(v);
end

% Masa de la densidad constante a trozos rho, definida sobre las celdas de x,
% acumulada en cada celda de la malla y.
function m = masa(x, rho, y)
  theta = [0, cumsum(diff(x).*rho)];
  t = interp1(x, theta, y, 'linear', 'extrap');
  m = diff(t);
end
