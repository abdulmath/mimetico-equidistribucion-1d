% prueba_dos_problemas.m
% Comprueba que el proceso de malla adaptada no está ajustado a un solo caso:
% se corre sobre dos problemas con capas de transición muy distintas y se mide
% convergencia del lazo, error final y orden en la sucesión de mallas.

function prueba_dos_problemas()
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

for P = [A B]
  printf('\n=== %s ===\n', P.nombre);

  printf('  Malla uniforme:\n');
  printf('    %-8s %-13s %-7s %-13s %-7s\n', 'Celdas', 'Linf', 'orden', 'L2', 'orden');
  eI0 = NaN; e20 = NaN;
  for N = [200 400 800 1600 3200]
    x = linspace(0, 1, N+1);
    U = resuelveRobin1D(x, P.k, P.f, P.ga, P.gb);
    [eI, e2] = erroresMalla(x, U, P.ue);
    imprime(N, eI, e2, eI0, e20);
    eI0 = eI;  e20 = e2;
  end

  printf('  Malla adaptada, mismo número de nodos:\n');
  printf('    %-8s %-6s %-6s %-13s %-7s %-13s %-7s %-9s\n', 'Celdas', 'iter', ...
         'conv', 'Linf', 'orden', 'L2', 'orden', 'ganancia');
  o = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);
  eI0 = NaN; e20 = NaN;
  for N = [200 400 800 1600 3200]
    xu = linspace(0, 1, N+1);
    [xa, Ua, h] = adaptaDeBoor(xu, P.k, P.f, P.ga, P.gb, o);
    [eI, e2] = erroresMalla(xa, Ua, P.ue);
    Uu = resuelveRobin1D(xu, P.k, P.f, P.ga, P.gb);
    [eIu, ~] = erroresMalla(xu, Uu, P.ue);
    if h.convergio, s = 'sí'; else, s = 'NO'; end
    if isnan(eI0)
      printf('    %-8d %-6d %-6s %-13.4e %-7s %-13.4e %-7s %-9.1f\n', ...
             N, h.it, s, eI, '--', e2, '--', eIu/eI);
    else
      printf('    %-8d %-6d %-6s %-13.4e %-7.2f %-13.4e %-7.2f %-9.1f\n', ...
             N, h.it, s, eI, log2(eI0/eI), e2, log2(e20/e2), eIu/eI);
    end
    eI0 = eI;  e20 = e2;
  end

  printf('  Efecto de la anchura de la ventana del filtro:\n');
  printf('    %-8s %-13s %-7s %-11s %-13s %-7s %-11s\n', 'Celdas', ...
         'Linf anchura', 'orden', 'salto', 'Linf D-D', 'orden', 'salto');
  oA = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);
  oD = oA;  oD.filtro = 'dorfidrury';
  e10 = NaN;  e30 = NaN;
  for N = [200 400 800 1600 3200]
    xu = linspace(0, 1, N+1);
    [xa, Ua] = adaptaDeBoor(xu, P.k, P.f, P.ga, P.gb, oA);
    [xd, Ud] = adaptaDeBoor(xu, P.k, P.f, P.ga, P.gb, oD);
    [e1, ~] = erroresMalla(xa, Ua, P.ue);
    [e3, ~] = erroresMalla(xd, Ud, P.ue);
    printf('    %-8d %-13.4e %-7s %-11.4f %-13.4e %-7s %-11.4f\n', N, ...
           e1, orden(e10, e1), salto(xa), e3, orden(e30, e3), salto(xd));
    e10 = e1;  e30 = e3;
  end
end
printf('\n');
end

% --------------------------------------------------------------------------
function s = orden(ant, act)
  if isnan(ant), s = '--'; else, s = sprintf('%.2f', log2(ant/act)); end
end

function s = salto(x)
  h = diff(x);
  s = max(abs(h(2:end)./h(1:end-1) - 1));
end

function imprime(N, eI, e2, eI0, e20)
  if isnan(eI0)
    printf('    %-8d %-13.4e %-7s %-13.4e %-7s\n', N, eI, '--', e2, '--');
  else
    printf('    %-8d %-13.4e %-7.2f %-13.4e %-7.2f\n', N, eI, log2(eI0/eI), e2, log2(e20/e2));
  end
end
