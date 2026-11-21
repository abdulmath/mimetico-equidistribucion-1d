% prueba_orden3.m
% Tercera parte: comprobar si el salto de malla que no decrece se debe a que
% la densidad de las dos celdas extremas se construye con los valores del
% gradiente que salen de la plantilla de frontera. Se compara el lazo tal cual
% con una variante que extrapola la densidad en esas dos celdas.

function prueba_orden3()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
kf = @(z) 1/al + al*(z-x0).^2;
ff = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
ga = al/(1+al^2*x0^2);
gb = -atan(al*(1-x0)) - atan(al*x0);

printf('\n=== Efecto de extrapolar la densidad en las celdas extremas ===\n');
printf('  %-7s | %-11s %-6s %-9s | %-11s %-6s %-9s\n', 'Celdas', ...
       'Linf tal cual', 'orden', 'salto', 'Linf extrap.', 'orden', 'salto');
p1 = NaN;  p2 = NaN;
for N = [100 200 400 800 1600]
  [e1, s1] = lazo(N, kf, ff, ga, gb, ue, false);
  [e2, s2] = lazo(N, kf, ff, ga, gb, ue, true);
  printf('  %-7d | %-11.4e %-6s %-9.4f | %-11.4e %-6s %-9.4f\n', N, ...
         e1, ord(p1,e1), s1, e2, ord(p2,e2), s2);
  p1 = e1;  p2 = e2;
end
printf('\n');
end

% --------------------------------------------------------------------------
function [eI, salto] = lazo(N, kf, ff, ga, gb, ue, extrap)
  x = linspace(0, 1, N+1);
  for j = 1:300
    [U, G] = resuelveRobin1D(x, kf, ff, ga, gb);
    gu  = G*U;
    dUc = (gu(1:N) + gu(2:N+1)).'/2;
    if extrap
      dUc(1)   = 2*dUc(2)     - dUc(3);
      dUc(end) = 2*dUc(end-1) - dUc(end-2);
    end
    rho = suavizaDorfiDrury(monitorArco(dUc, 0.5, log(1.052), 1), 3, 2);
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

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
