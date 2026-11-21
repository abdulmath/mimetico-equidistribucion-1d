% prueba_orden7.m
% Séptima parte: la anchura local sigma_i = max(sigma0, C h_i) fue la única
% variante que convergió en casi toda la sucesión, y falló sólo en la malla
% más gruesa del problema A. Aquí se barre C para ver si existe un valor que
% haga converger toda la sucesión en los dos problemas sin perder el orden.

function prueba_orden7()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
A.nombre = 'A';
A.k  = @(z) 1/al + al*(z-x0).^2;
A.f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
A.ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
A.ga = al/(1+al^2*x0^2);
A.gb = -atan(al*(1-x0)) - atan(al*x0);

c = atan(100);
B.nombre = 'B';
B.k  = @(z) ones(size(z));
B.f  = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue = @(z) atan(100*z)/c;
B.ga = 100/c;
B.gb = 1 + 100/(c*(1+1e4));

Cs = [2 3 4 6 8];
printf('\n=== Anchura local sigma_i = max(0.02, C h_i): barrido de C ===\n');
for P = [A B]
  for ci = 1:numel(Cs)
    printf('  Problema %s, C = %g\n', P.nombre, Cs(ci));
    printf('    %-7s %-6s %-13s %-7s %-11s %-11s\n', 'Celdas', 'iter', ...
           'Linf', 'orden', 'salto', 'h_max/h_min');
    p = NaN;
    for N = [100 200 400 800 1600]
      r = lazo(N, P, 0.02, Cs(ci));
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
function r = lazo(N, P, sigma0, C)
  x = linspace(0, 1, N+1);
  r.convergio = false;  r.it = 300;
  for j = 1:300
    [U, G] = resuelveRobin1D(x, P.k, P.f, P.ga, P.gb);
    gu  = G*U;
    dUc = (gu(1:N) + gu(2:N+1)).'/2;
    rho = monitorArco(dUc, 0.5, log(1.052), 1);
    rho = suavizaLocal(x, rho, sigma0, C);
    xn  = deBoor(x, rho);
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

function rs = suavizaLocal(x, rho, sigma0, C)
  rho = rho(:).';
  xc  = (x(1:end-1) + x(2:end))/2;
  h   = diff(x);
  sig = max(sigma0, C*h(:));
  d   = (xc.' - xc) ./ sig;
  w   = exp(-d.^2) .* h;
  rs  = ((w*rho.') ./ sum(w, 2)).';
end

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
