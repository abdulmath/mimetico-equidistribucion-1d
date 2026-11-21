% prueba_tolerancia.m
% Por qué la comparación del suavizado elíptico se detiene con tolX = 1e-9 y
% no con el 1e-10 que usa el resto del artículo.
%
% Con 1e-10 el problema A con 3200 celdas y beta_s = 20 no converge: la
% iteración exterior entra en un ciclo y el movimiento de la malla oscila sin
% volver a descender. El error final no cambia en las cifras que se reportan,
% de modo que los valores de la tabla son los mismos por cualquiera de las dos
% vías, pero la malla del suavizado elíptico queda determinada un orden de
% magnitud peor que la del filtro de anchura fija, que sí alcanza 1e-10.
%
% Tarda unos minutos: la corrida que no converge agota las 300 iteraciones.

function prueba_tolerancia()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
k  = @(z) 1/al + al*(z-x0).^2;
f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
ga = al/(1+al^2*x0^2);
gb = -atan(al*(1-x0)) - atan(al*x0);

N  = 3200;
xu = linspace(0, 1, N+1);

V = { 'global, betaS=20', struct('filtro','global','betaS',20)
      'anchura fija',     struct('filtro','anchura')          };

printf('\n=== Problema A, %d celdas, dos tolerancias ===\n', N);
printf('  %-20s %-9s %-7s %-7s %-13s %-13s\n', 'suavizado', 'tolX', ...
       'conv', 'iter', 'mov. final', 'Linf');

hist10 = [];
for v = 1:size(V,1)
  for tol = [1e-9 1e-10]
    o = V{v,2};
    o.omega = 1;  o.maxIt = 300;  o.tolX = tol;
    [xa, Ua, h] = adaptaDeBoor(xu, k, f, ga, gb, o);
    eI = erroresMalla(xa, Ua, ue);
    if h.convergio, s = 'sí'; else, s = 'no'; end
    printf('  %-20s %-9.0e %-7s %-7d %-13.3e %-13.4e\n', ...
           V{v,1}, tol, s, h.it, h.mov(end), eI);
    if v == 1 && tol == 1e-10, hist10 = h.mov; end
    fflush(stdout);
  end
end

if ~isempty(hist10) && numel(hist10) > 60
  cola = hist10(51:end);
  printf('\n=== El ciclo del suavizado elíptico con tolX = 1e-10 ===\n');
  printf('  Movimiento de la malla a partir de la iteración 50.\n');
  printf('  mínimo  %.3e\n', min(cola));
  printf('  máximo  %.3e\n', max(cola));
  printf('  mediana %.3e\n', median(cola));
  printf('\n  Últimas diez iteraciones\n');
  for j = numel(hist10)-9:numel(hist10)
    printf('  %-6d %.3e\n', j, hist10(j));
  end
end
printf('\n');
end
