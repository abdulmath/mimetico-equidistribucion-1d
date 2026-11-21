% prueba_orden6.m
% Sexta parte: arreglar el único defecto que dejó el filtro de anchura fija en
% x, a saber, que en mallas gruesas el punto fijo no converge. Se prueban tres
% remedios: hacer la anchura del núcleo adaptativa de forma global, hacerla
% adaptativa celda a celda, y amortiguar el lazo con relajación. Se busca la
% variante que converja en toda la sucesión de mallas sin perder el orden 2.

function prueba_orden6()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
A.nombre = 'A: pico interior en x0 = 0.8';
A.k  = @(z) 1/al + al*(z-x0).^2;
A.f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
A.ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
A.ga = al/(1+al^2*x0^2);
A.gb = -atan(al*(1-x0)) - atan(al*x0);

c = atan(100);
B.nombre = 'B: capa en x = 0';
B.k  = @(z) ones(size(z));
B.f  = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue = @(z) atan(100*z)/c;
B.ga = 100/c;
B.gb = 1 + 100/(c*(1+1e4));

% Cada variante: {etiqueta, modo, sigma0, C, omega}
%   modo 'fijo'   anchura constante sigma0
%   modo 'global' anchura max(sigma0, C*max(h)), una sola para toda la malla
%   modo 'local'  anchura max(sigma0, C*h_i), distinta en cada celda
V = { 'fijo 0.02',        'fijo',   0.02, 0, 1.0
      'global C=3',       'global', 0.02, 3, 1.0
      'local C=3',        'local',  0.02, 3, 1.0
      'fijo 0.02 w=0.5',  'fijo',   0.02, 0, 0.5 };

for P = [A B]
  printf('\n=== Problema %s ===\n', P.nombre);
  for v = 1:size(V,1)
    printf('  %s\n', V{v,1});
    printf('    %-7s %-6s %-13s %-7s %-11s %-11s\n', 'Celdas', 'iter', ...
           'Linf', 'orden', 'salto', 'h_max/h_min');
    p = NaN;
    for N = [100 200 400 800 1600]
      r = lazo(N, P, 0.5, V{v,2}, V{v,3}, V{v,4}, V{v,5});
      if r.convergio, et = sprintf('%d', r.it); else, et = sprintf('%d!', r.it); end
      printf('    %-7d %-6s %-13.4e %-7s %-11.4f %-11.1f\n', ...
             N, et, r.eI, ord(p, r.eI), r.salto, r.razon);
      p = r.eI;
    end
  end
end
printf('\n  El signo ! marca que el lazo agotó las iteraciones sin converger.\n\n');
end

% --------------------------------------------------------------------------
function r = lazo(N, P, beta, modo, sigma0, C, omega)
  x = linspace(0, 1, N+1);
  r.convergio = false;  r.it = 300;
  for j = 1:300
    [U, G] = resuelveRobin1D(x, P.k, P.f, P.ga, P.gb);
    gu  = G*U;
    dUc = (gu(1:N) + gu(2:N+1)).'/2;
    rho = monitorArco(dUc, beta, log(1.052), 1);
    rho = suavizaAnchuraFija(x, rho, modo, sigma0, C);
    xeq = deBoor(x, rho);
    xn  = (1-omega)*x + omega*xeq;
    xn(1) = 0;  xn(end) = 1;
    mov = max(abs(xn - x));
    x   = xn;
    if mov < 1e-10, r.convergio = true;  r.it = j;  break, end
  end
  U = resuelveRobin1D(x, P.k, P.f, P.ga, P.gb);
  [r.eI, r.e2] = erroresMalla(x, U, P.ue);
  h = diff(x);
  r.salto = max(abs(h(2:end)./h(1:end-1) - 1));
  r.razon = max(h)/min(h);
end

% Promedio con núcleo gaussiano en la variable x, pesado por la longitud de
% cada celda. La anchura puede ser constante, global adaptativa o local.
function rs = suavizaAnchuraFija(x, rho, modo, sigma0, C)
  rho = rho(:).';
  xc  = (x(1:end-1) + x(2:end))/2;
  h   = diff(x);
  switch modo
    case 'fijo',   sig = sigma0;
    case 'global', sig = max(sigma0, C*max(h));
    case 'local',  sig = max(sigma0, C*h(:));
  end
  d  = (xc.' - xc) ./ sig;
  w  = exp(-d.^2) .* h;
  rs = ((w*rho.') ./ sum(w, 2)).';
end

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
