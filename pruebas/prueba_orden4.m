% prueba_orden4.m
% Cuarta parte: contrastar la explicación que dejaron las pruebas anteriores.
% Con densidad fija e independiente de la malla el orden vuelve a 2, así que
% lo que estropea el orden es que el monitor se construya sobre la propia
% malla en movimiento. El filtro de Dorfi y Drury promedia un número fijo de
% celdas, de modo que su ventana se encoge con la malla y nunca llega a
% suavizar la densidad como función de x. Aquí se sustituye ese filtro por uno
% de anchura fija en x y se mira si el salto de malla vuelve a decrecer.

function prueba_orden4()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
kf = @(z) 1/al + al*(z-x0).^2;
ff = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
ga = al/(1+al^2*x0^2);
gb = -atan(al*(1-x0)) - atan(al*x0);

printf('\n=== Filtro de anchura fija en x frente al de anchura fija en celdas ===\n');
printf('  %-7s | %-11s %-6s %-8s | %-11s %-6s %-8s | %-11s %-6s %-8s\n', ...
       'Celdas', 'Dorfi-Drury', 'orden', 'salto', 'sigma=0.02', 'orden', ...
       'salto', 'sigma=0.05', 'orden', 'salto');
p1 = NaN;  p2 = NaN;  p3 = NaN;
for N = [100 200 400 800 1600]
  [e1, s1] = lazo(N, kf, ff, ga, gb, ue, 0.5, 0);
  [e2, s2] = lazo(N, kf, ff, ga, gb, ue, 0.5, 0.02);
  [e3, s3] = lazo(N, kf, ff, ga, gb, ue, 0.5, 0.05);
  printf('  %-7d | %-11.4e %-6s %-8.4f | %-11.4e %-6s %-8.4f | %-11.4e %-6s %-8.4f\n', ...
         N, e1, ord(p1,e1), s1, e2, ord(p2,e2), s2, e3, ord(p3,e3), s3);
  p1 = e1;  p2 = e2;  p3 = e3;
end

printf('\n=== El mismo filtro fijo, con menos escala en el monitor ===\n');
printf('  %-7s | %-11s %-6s %-8s | %-11s %-6s %-8s\n', 'Celdas', ...
       'beta=0.1 s=0.02', 'orden', 'salto', 'beta=0.1 DD', 'orden', 'salto');
q1 = NaN;  q2 = NaN;
for N = [100 200 400 800 1600]
  [e1, s1] = lazo(N, kf, ff, ga, gb, ue, 0.1, 0.02);
  [e2, s2] = lazo(N, kf, ff, ga, gb, ue, 0.1, 0);
  printf('  %-7d | %-11.4e %-6s %-8.4f | %-11.4e %-6s %-8.4f\n', ...
         N, e1, ord(q1,e1), s1, e2, ord(q2,e2), s2);
  q1 = e1;  q2 = e2;
end
printf('\n');
end

% --------------------------------------------------------------------------
% Lazo de malla en movimiento. sigma = 0 usa el filtro de Dorfi y Drury;
% sigma > 0 usa un promedio gaussiano de anchura fija en x.
function [eI, salto] = lazo(N, kf, ff, ga, gb, ue, beta, sigma)
  x = linspace(0, 1, N+1);
  for j = 1:300
    [U, G] = resuelveRobin1D(x, kf, ff, ga, gb);
    gu  = G*U;
    dUc = (gu(1:N) + gu(2:N+1)).'/2;
    rho = monitorArco(dUc, beta, log(1.052), 1);
    if sigma > 0
      rho = suavizaAnchuraFija(x, rho, sigma);
    else
      rho = suavizaDorfiDrury(rho, 3, 2);
    end
    xn  = deBoor(x, rho);
    mov = max(abs(xn - x));
    x   = xn;
    if mov < 1e-10, break, end
  end
  U = resuelveRobin1D(x, kf, ff, ga, gb);
  [eI, ~] = erroresMalla(x, U, ue);
  h = diff(x);
  salto = max(abs(h(2:end)./h(1:end-1) - 1));
end

% Promedio con núcleo gaussiano de anchura sigma en la variable x, pesado por
% la longitud de cada celda para que no dependa de cómo estén repartidas.
function rs = suavizaAnchuraFija(x, rho, sigma)
  rho = rho(:).';
  xc  = (x(1:end-1) + x(2:end))/2;
  h   = diff(x);
  d   = (xc.' - xc)/sigma;
  w   = exp(-d.^2) .* h;
  rs  = (w*rho.') ./ sum(w, 2);
  rs  = rs.';
end

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
