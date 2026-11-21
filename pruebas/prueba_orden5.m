% prueba_orden5.m
% Quinta parte: comprobar que la conclusión de prueba_orden4 no depende de un
% valor afortunado de la anchura del filtro ni del problema concreto. Se barre
% sigma sobre el problema A y se repite el problema B, que ya daba orden 2 con
% el filtro de Dorfi y Drury, para confirmar que el filtro de anchura fija no
% lo estropea.

function prueba_orden5()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

al = 250;  x0 = 0.8;
A.k  = @(z) 1/al + al*(z-x0).^2;
A.f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
A.ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
A.ga = al/(1+al^2*x0^2);
A.gb = -atan(al*(1-x0)) - atan(al*x0);

c = atan(100);
B.k  = @(z) ones(size(z));
B.f  = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue = @(z) atan(100*z)/c;
B.ga = 100/c;
B.gb = 1 + 100/(c*(1+1e4));

printf('\n=== Problema A: barrido de la anchura del filtro ===\n');
printf('  La capa del problema A tiene anchura del orden de 1/a = %.4f\n', 1/al);
printf('  %-7s | %-10s %-6s | %-10s %-6s | %-10s %-6s | %-10s %-6s\n', 'Celdas', ...
       'sigma=0.01', 'orden', 'sigma=0.02', 'orden', 'sigma=0.05', 'orden', ...
       'sigma=0.10', 'orden');
p = NaN(1,4);
for N = [100 200 400 800 1600]
  e = zeros(1,4);  s = [0.01 0.02 0.05 0.10];
  for m = 1:4, [e(m), ~] = lazo(N, A, 0.5, s(m)); end
  printf('  %-7d | %-10.4e %-6s | %-10.4e %-6s | %-10.4e %-6s | %-10.4e %-6s\n', ...
         N, e(1), ord(p(1),e(1)), e(2), ord(p(2),e(2)), e(3), ord(p(3),e(3)), ...
         e(4), ord(p(4),e(4)));
  p = e;
end

printf('\n=== Problema B: el filtro de anchura fija no estropea lo que ya iba bien ===\n');
printf('  %-7s | %-11s %-6s %-8s | %-11s %-6s %-8s\n', 'Celdas', ...
       'Dorfi-Drury', 'orden', 'salto', 'sigma=0.02', 'orden', 'salto');
q1 = NaN;  q2 = NaN;
for N = [100 200 400 800 1600]
  [e1, s1] = lazo(N, B, 0.5, 0);
  [e2, s2] = lazo(N, B, 0.5, 0.02);
  printf('  %-7d | %-11.4e %-6s %-8.4f | %-11.4e %-6s %-8.4f\n', ...
         N, e1, ord(q1,e1), s1, e2, ord(q2,e2), s2);
  q1 = e1;  q2 = e2;
end

printf('\n=== Iteraciones y razón de tamaños con el filtro de anchura fija ===\n');
printf('  %-9s %-7s | %-7s %-11s | %-7s %-11s\n', 'Problema', 'Celdas', ...
       'iter', 'h_max/h_min', 'iter', 'h_max/h_min');
printf('  %-9s %-7s | %-19s | %-19s\n', '', '', 'Dorfi-Drury', 'sigma=0.02');
for N = [200 800 1600]
  [~, ~, i1, r1] = lazo(N, A, 0.5, 0);
  [~, ~, i2, r2] = lazo(N, A, 0.5, 0.02);
  printf('  %-9s %-7d | %-7d %-11.1f | %-7d %-11.1f\n', 'A', N, i1, r1, i2, r2);
end
for N = [200 800 1600]
  [~, ~, i1, r1] = lazo(N, B, 0.5, 0);
  [~, ~, i2, r2] = lazo(N, B, 0.5, 0.02);
  printf('  %-9s %-7d | %-7d %-11.1f | %-7d %-11.1f\n', 'B', N, i1, r1, i2, r2);
end
printf('\n');
end

% --------------------------------------------------------------------------
function [eI, salto, it, razon] = lazo(N, P, beta, sigma)
  x = linspace(0, 1, N+1);
  it = 300;
  for j = 1:300
    [U, G] = resuelveRobin1D(x, P.k, P.f, P.ga, P.gb);
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
    if mov < 1e-10, it = j; break, end
  end
  U = resuelveRobin1D(x, P.k, P.f, P.ga, P.gb);
  [eI, ~] = erroresMalla(x, U, P.ue);
  h = diff(x);
  salto = max(abs(h(2:end)./h(1:end-1) - 1));
  razon = max(h)/min(h);
end

function rs = suavizaAnchuraFija(x, rho, sigma)
  rho = rho(:).';
  xc  = (x(1:end-1) + x(2:end))/2;
  h   = diff(x);
  d   = (xc.' - xc)/sigma;
  w   = exp(-d.^2) .* h;
  rs  = ((w*rho.') ./ sum(w, 2)).';
end

function s = ord(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end
