% prueba_orden.m
% Indaga por qué el orden sobre mallas adaptadas se queda por debajo de 2 en
% el problema del pico interior. Se mira dónde está el error máximo, cómo de
% suave es la malla, y qué pasa al quitar la transformación de Tapia y López.

function prueba_orden()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
kf = @(z) 1/al + al*(z-x0).^2;
ff = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
ga = al/(1+al^2*x0^2);
gb = -atan(al*(1-x0)) - atan(al*x0);

% --- 1. dónde está el error máximo ---------------------------------------
printf('\n=== Dónde cae el error máximo, malla adaptada ===\n');
printf('  %-7s %-12s %-9s %-12s %-12s %-12s\n', 'Celdas', 'Linf', 'posición', ...
       'Linf interior', 'h_1/h_med', 'h_N/h_med');
o = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);
for N = [100 200 400 800 1600]
  xu = linspace(0, 1, N+1);
  [xa, Ua] = adaptaDeBoor(xu, kf, ff, ga, gb, o);
  [e, u, h] = residuo(xa, Ua, ue);
  [eI, k] = max(abs(e));
  eInt = max(abs(e(3:end-2)));
  printf('  %-7d %-12.4e %-9.4f %-12.4e %-12.3f %-12.3f\n', ...
         N, eI, u(k), eInt, h(1)/mean(h), h(end)/mean(h));
end

% --- 2. suavidad de la malla ---------------------------------------------
printf('\n=== Suavidad de la malla adaptada ===\n');
printf('  Se mide max|h_{i+1}/h_i - 1|; la supraconvergencia del esquema pide\n');
printf('  que esa cantidad tienda a cero al refinar.\n');
printf('  %-7s %-14s %-14s %-14s\n', 'Celdas', 'max salto', 'max salto s/TL', 'h_max/h_min');
for N = [100 200 400 800 1600]
  xu = linspace(0, 1, N+1);
  [xa, ~] = adaptaDeBoor(xu, kf, ff, ga, gb, o);
  h  = diff(xa);
  o2 = o;  o2.k1 = [];  o2.k2 = 0;
  [xb, ~] = adaptaDeBoor(xu, kf, ff, ga, gb, o2);
  hb = diff(xb);
  printf('  %-7d %-14.4f %-14.4f %-14.1f\n', N, ...
         max(abs(h(2:end)./h(1:end-1) - 1)), ...
         max(abs(hb(2:end)./hb(1:end-1) - 1)), max(h)/min(h));
end

% --- 3. orden con y sin la transformación de Tapia y López ---------------
printf('\n=== Orden según la densidad empleada ===\n');
printf('  %-7s %-12s %-7s %-12s %-7s %-12s %-7s\n', 'Celdas', ...
       'Linf con TL', 'orden', 'Linf sin TL', 'orden', 'Linf arco b=5', 'orden');
p1 = NaN; p2 = NaN; p3 = NaN;
oA = o;                        % con Tapia-López
oB = o;  oB.k1 = [];  oB.k2 = 0;             % sólo longitud de arco
oC = oB; oC.beta = 5;                        % arco con más escala
for N = [100 200 400 800 1600]
  xu = linspace(0, 1, N+1);
  e1 = errAdapt(xu, kf, ff, ga, gb, oA, ue);
  e2 = errAdapt(xu, kf, ff, ga, gb, oB, ue);
  e3 = errAdapt(xu, kf, ff, ga, gb, oC, ue);
  printf('  %-7d %-12.4e %-7s %-12.4e %-7s %-12.4e %-7s\n', N, ...
         e1, ord(p1,e1), e2, ord(p2,e2), e3, ord(p3,e3));
  p1 = e1;  p2 = e2;  p3 = e3;
end
printf('\n');
end

% --------------------------------------------------------------------------
function [e, u, h] = residuo(x, U, ue)
  N = numel(x) - 1;
  u = [x(1), (x(1:N) + x(2:N+1))/2, x(end)];
  e = ue(u) - U(:).';
  h = diff(x);
end

function eI = errAdapt(xu, kf, ff, ga, gb, o, ue)
  [xa, Ua] = adaptaDeBoor(xu, kf, ff, ga, gb, o);
  [eI, ~] = erroresMalla(xa, Ua, ue);
end

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
