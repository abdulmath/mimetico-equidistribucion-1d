% prueba_orden2.m
% Segunda parte de la indagación: localizar el salto de malla que no decrece
% y separar dos causas posibles del orden bajo, a saber, que la densidad se
% calcule a partir del gradiente discreto sobre la propia malla que se está
% moviendo, o que la densidad de longitud de arco no sea la adecuada.

function prueba_orden2()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
kf = @(z) 1/al + al*(z-x0).^2;
ff = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
du = @(z) -(atan(al*(z-x0)) + atan(al*x0)) + (1-z).*(al./(1+al^2*(z-x0).^2));
ga = al/(1+al^2*x0^2);
gb = -atan(al*(1-x0)) - atan(al*x0);
o  = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);

% --- 1. dónde está el salto que no decrece --------------------------------
printf('\n=== Localización del mayor salto de malla ===\n');
printf('  %-7s %-9s %-9s %-11s %-11s %-11s\n', 'Celdas', 'índice', 'posición', ...
       'salto max', 'salto medio', 'salto interior');
for N = [200 400 800 1600]
  xu = linspace(0, 1, N+1);
  [xa, ~] = adaptaDeBoor(xu, kf, ff, ga, gb, o);
  h = diff(xa);
  r = abs(h(2:end)./h(1:end-1) - 1);
  [rm, k] = max(r);
  printf('  %-7d %-9d %-9.4f %-11.4f %-11.4f %-11.4f\n', N, k, xa(k+1), rm, ...
         mean(r), max(r(4:end-3)));
end

% --- 2. malla generada por una densidad fija e independiente de la malla ---
printf('\n=== Malla equidistribuida con densidad fija (derivada exacta) ===\n');
printf('  La densidad ya no depende de la malla: si el orden sube a 2, el\n');
printf('  responsable es el monitor calculado sobre la malla en movimiento.\n');
printf('  %-7s %-13s %-7s %-13s %-7s %-11s\n', 'Celdas', 'Linf', 'orden', ...
       'L2', 'orden', 'salto max');
p1 = NaN;  p2 = NaN;
for N = [100 200 400 800 1600]
  x = mallaFija(N, du, 0.5, 4000);
  U = resuelveRobin1D(x, kf, ff, ga, gb);
  [eI, e2] = erroresMalla(x, U, ue);
  h = diff(x);
  printf('  %-7d %-13.4e %-7s %-13.4e %-7s %-11.4f\n', N, eI, ord(p1,eI), ...
         e2, ord(p2,e2), max(abs(h(2:end)./h(1:end-1) - 1)));
  p1 = eI;  p2 = e2;
end

% --- 3. la misma densidad fija, sin suavizar y con varias escalas ---------
printf('\n=== Densidad fija: efecto del parámetro de escala ===\n');
printf('  %-7s %-13s %-7s %-13s %-7s %-13s %-7s\n', 'Celdas', ...
       'beta = 0.1', 'orden', 'beta = 0.5', 'orden', 'beta = 5', 'orden');
q1 = NaN;  q2 = NaN;  q3 = NaN;
for N = [100 200 400 800 1600]
  e1 = errFija(N, du, 0.1, kf, ff, ga, gb, ue);
  e2 = errFija(N, du, 0.5, kf, ff, ga, gb, ue);
  e3 = errFija(N, du, 5.0, kf, ff, ga, gb, ue);
  printf('  %-7d %-13.4e %-7s %-13.4e %-7s %-13.4e %-7s\n', N, ...
         e1, ord(q1,e1), e2, ord(q2,e2), e3, ord(q3,e3));
  q1 = e1;  q2 = e2;  q3 = e3;
end
printf('\n');
end

% --------------------------------------------------------------------------
% Malla que equidistribuye la longitud de arco de la derivada exacta,
% calculada sobre una malla auxiliar muy fina para que no dependa de N.
function x = mallaFija(N, du, beta, M)
  xf  = linspace(0, 1, M+1);
  xc  = (xf(1:M) + xf(2:M+1))/2;
  rho = sqrt(1 + beta*du(xc).^2);
  x   = deBoor(xf, rho, N);
end

function eI = errFija(N, du, beta, kf, ff, ga, gb, ue)
  x = mallaFija(N, du, beta, 4000);
  U = resuelveRobin1D(x, kf, ff, ga, gb);
  [eI, ~] = erroresMalla(x, U, ue);
end

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
