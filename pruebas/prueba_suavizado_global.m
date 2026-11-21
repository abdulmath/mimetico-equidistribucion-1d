% prueba_suavizado_global.m
% Compara, dentro del mismo lazo de de Boor y sobre los dos problemas de
% prueba, tres formas de suavizar la función monitor:
%
%   'anchura'     filtro gaussiano de anchura fija en la variable espacial x
%   'dorfidrury'  filtro de pesos geométricos sobre un número fijo de celdas
%   'global'      suavizado elíptico (I - beta^{-2} d^2/dxi^2) rs = rho
%                 resuelto en la variable computacional xi, con beta fijo
%
% La tercera opción es la de Huang y Russell, de la que el filtro de Dorfi y
% Drury es la truncación local. Su teorema para el esquema en diferencias
% acota
%
%   nu <= rs_j/rs_{j-1} <= 1/nu,   nu = (sqrt(1+4/(beta^2 dxi^2)) - 1)
%                                     /(sqrt(1+4/(beta^2 dxi^2)) + 1),
%
% y con beta fijo dxi = 1/N tiende a cero, luego nu tiende a uno. La columna
% nu de la salida permite contrastar esa cota contra el salto medido.
%
% Se imprimen error, orden, max|h_{i+1}/h_i - 1| y número de iteraciones
% exteriores. Al final se mide el coste por llamada de los tres filtros.

function prueba_suavizado_global()
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

Ns    = [200 400 800 1600 3200];
betas = [5 10 20 50];

for P = [A B]
  printf('\n================ %s ================\n', P.nombre);

  printf('\n-- anchura fija en x (sigma0 = 0.02, C = 2) --\n');
  corre(P, Ns, struct('omega',1,'maxIt',300,'tolX',1e-9,'filtro','anchura'));

  printf('\n-- Dorfi y Drury (p = 3, q = 2) --\n');
  corre(P, Ns, struct('omega',1,'maxIt',300,'tolX',1e-9,'filtro','dorfidrury'));

  for bs = betas
    printf('\n-- global elíptico en xi, beta = %g --\n', bs);
    corre(P, Ns, struct('omega',1,'maxIt',300,'tolX',1e-9, ...
                        'filtro','global','betaS',bs));
  end
end

printf('\n================ coste por llamada al filtro ================\n');
printf('   %-8s %-12s %-12s %-12s\n', 'N', 'anchura', 'global', 'Dorfi-Drury');
for N = [400 800 1600 3200 6400]
  x   = linspace(0, 1, N+1);
  rho = 1 + 0.5*sin(10*pi*(x(1:N) + x(2:N+1))/2);
  t = tic;  for r = 1:3, suavizaAnchuraFija(x, rho, 0.02, 2); end;  ta = toc(t)/3;
  t = tic;  for r = 1:3, suavizaGlobal(rho, 20);              end;  tg = toc(t)/3;
  t = tic;  for r = 1:3, suavizaDorfiDrury(rho, 3, 2);        end;  td = toc(t)/3;
  printf('   %-8d %-12.5f %-12.5f %-12.5f\n', N, ta, tg, td);
  fflush(stdout);
end
printf('\n');
end

% --------------------------------------------------------------------------
function corre(P, Ns, o)
  printf('   %-7s %-6s %-5s %-12s %-6s %-12s %-6s %-10s %-8s\n', ...
         'N', 'iter', 'conv', 'Linf', 'orden', 'L2', 'orden', 'salto', 'nu');
  e10 = NaN;  e20 = NaN;
  for N = Ns
    xu = linspace(0, 1, N+1);
    [xa, Ua, h] = adaptaDeBoor(xu, P.k, P.f, P.ga, P.gb, o);
    [e1, e2] = erroresMalla(xa, Ua, P.ue);
    if h.convergio, s = 'sí'; else, s = 'NO'; end
    if strcmp(o.filtro, 'global')
      dxi = 1/N;  t = sqrt(1 + 4/(o.betaS^2*dxi^2));
      nu  = sprintf('%.4f', (t-1)/(t+1));
    else
      nu = '--';
    end
    if isnan(e10)
      o1 = '--';  o2 = '--';
    else
      o1 = sprintf('%.2f', log2(e10/e1));  o2 = sprintf('%.2f', log2(e20/e2));
    end
    hh = diff(xa);  sa = max(abs(hh(2:end)./hh(1:end-1) - 1));
    printf('   %-7d %-6d %-5s %-12.4e %-6s %-12.4e %-6s %-10.4f %-8s\n', ...
           N, h.it, s, e1, o1, e2, o2, sa, nu);
    fflush(stdout);
    e10 = e1;  e20 = e2;
  end
end
