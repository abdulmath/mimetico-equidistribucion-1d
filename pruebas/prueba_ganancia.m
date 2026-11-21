% prueba_ganancia.m
% Prolonga la sucesión de mallas una etapa más allá de las tablas del
% artículo, hasta 6400 celdas, para comprobar que la ganancia frente a la
% malla uniforme se estabiliza en lugar de degradarse al refinar.
%
% Si el límite de la densidad suavizada fuese constante, la malla límite se
% aplanaría y la ganancia caería hacia uno. Lo que se observa es lo contrario:
% con 6400 celdas sigue valiendo lo mismo que con 3200.
%
% Tarda unos minutos. El filtro de escala fija cuesta O(N^2) y con 6400 celdas
% una sola llamada ronda el segundo.

function prueba_ganancia()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

% --- Problema A : pico interior, coeficiente variable ---------------------
al = 250;  x0 = 0.8;
A.nombre = 'A: pico interior en x0 = 0.8, k(x) = 1/a + a(x-x0)^2';
A.k  = @(z) 1/al + al*(z-x0).^2;
A.f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
A.ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
A.ga = al/(1+al^2*x0^2);
A.gb = -atan(al*(1-x0)) - atan(al*x0);

% --- Problema B : capa en el borde izquierdo, coeficiente constante -------
c = atan(100);
B.nombre = 'B: capa en x = 0, k(x) = 1';
B.k  = @(z) ones(size(z));
B.f  = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue = @(z) atan(100*z)/c;
B.ga = 100/c;
B.gb = 1 + 100/(c*(1+1e4));

o  = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);
Ns = [3200 6400];

for P = [A B]
  printf('\n=== %s ===\n', P.nombre);
  printf('  %-8s %-6s %-6s %-13s %-7s %-13s %-7s %-9s\n', 'Celdas', ...
         'iter', 'conv', 'Linf', 'orden', 'L2', 'orden', 'ganancia');
  eIant = NaN;  e2ant = NaN;
  for N = Ns
    xu = linspace(0, 1, N+1);
    [xa, Ua, h] = adaptaDeBoor(xu, P.k, P.f, P.ga, P.gb, o);
    [eIa, e2a]  = erroresMalla(xa, Ua, P.ue);

    Uu          = resuelveRobin1D(xu, P.k, P.f, P.ga, P.gb);
    [eIu, ~]    = erroresMalla(xu, Uu, P.ue);

    if h.convergio, s = 'sí'; else, s = 'no'; end
    printf('  %-8d %-6d %-6s %-13.4e %-7s %-13.4e %-7s %-9.1f\n', ...
           N, h.it, s, eIa, ord(eIant, eIa), e2a, ord(e2ant, e2a), eIu/eIa);
    eIant = eIa;  e2ant = e2a;
    fflush(stdout);
  end
end
printf('\n  La ganancia con 6400 celdas debe coincidir con la de 3200.\n\n');
end

% --------------------------------------------------------------------------
function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
