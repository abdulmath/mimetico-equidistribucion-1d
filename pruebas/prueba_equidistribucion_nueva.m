% prueba_equidistribucion_nueva.m
% Comprueba que el proceso de malla adaptada en movimiento converge: se
% imprime el movimiento máximo de la malla y el error en cada iteración
% exterior, y luego una tabla de errores finales para varios tamaños de malla.

function prueba_equidistribucion_nueva()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

a = 0;  b = 1;  alpha = 250;  x0 = 0.8;
kf = @(z) 1/alpha + alpha*(z-x0).^2;
ff = @(z) -2*(1 + alpha*(z-x0).*(atan(alpha*x0) + atan(alpha*(z-x0))));
ue = @(z) (1-z).*(atan(alpha*(z-x0)) + atan(alpha*x0));
ga = alpha/(1+alpha^2*x0^2);
gb = -atan(alpha*(1-x0)) - atan(alpha*x0);

opts = struct('beta', 0.5, 'k2', 1, 'nCont', 10, 'mInt', 3, 'maxIt', 60);

printf('\n=== Historia de la iteración exterior, 100 celdas ===\n');
printf('  %-5s %-13s %-13s %-13s\n', 'Iter', 'mov. malla', 'Linf', 'L2');
nelem = 100;  h = (b-a)/nelem;  x = a:h:b;
for j = 1:30
  o = opts;  o.maxIt = 1;  o.nCont = opts.nCont;
  % se reproduce la rampa a mano para poder imprimir cada paso
  [U, G] = resuelveRobin1D(x, kf, ff, ga, gb);
  [eI, e2] = erroresMalla(x, U, ue);
  gu  = G*U;
  dUc = (gu(1:nelem) + gu(2:nelem+1)).'/2;
  t   = min(1, j/opts.nCont);
  rho = monitorArco(dUc, opts.beta*t, log(1.052), opts.k2*t);
  xn  = equidistribuye(x, rho, opts.mInt, 3, 2);
  printf('  %-5d %-13.4e %-13.4e %-13.4e\n', j, max(abs(xn-x)), eI, e2);
  x = xn;
end

printf('\n=== Errores en la malla adaptada final ===\n');
printf('  %-8s %-8s %-8s %-13s %-13s %-13s %-13s\n', 'Celdas', 'nodos', ...
       'Iter', 'Linf adapt', 'L2 adapt', 'Linf unif', 'L2 unif');
for nelem = [40 60 80 100 110 120 160 200]
  h = (b-a)/nelem;  x0m = a:h:b;
  [xa, Ua, hist] = adaptaEquidistribucion(x0m, kf, ff, ga, gb, opts);
  [eIa, e2a] = erroresMalla(xa, Ua, ue);
  Uu = resuelveRobin1D(x0m, kf, ff, ga, gb);
  [eIu, e2u] = erroresMalla(x0m, Uu, ue);
  printf('  %-8d %-8d %-8d %-13.4e %-13.4e %-13.4e %-13.4e\n', ...
         nelem, nelem+1, hist.it, eIa, e2a, eIu, e2u);
end
printf('\n');
end
