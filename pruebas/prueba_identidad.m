% prueba_identidad.m
% Identidad mimética de los operadores del artículo sobre mallas no uniformes.
%
% Con Q = diag(h_1, h_1, ..., h_N, h_N) los volúmenes de la malla mimética, se
% busca la matriz P sobre las caras que hace exacta
%
%     D' Q + P G = B'
%
% Se resuelve P = (B' - D'Q)/G y se comprueba tres cosas: que P sale diagonal,
% que sale positiva, y qué relación guarda con los pesos de Newton-Cotes
% (3/8, 9/8, 1, ..., 1, 9/8, 3/8) escalados por la medida de cada cara,
% hf = (h_1, (h_1+h_2)/2, ..., (h_{N-1}+h_N)/2, h_N).

function prueba_identidad()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

printf('\n=== Identidad  D'' Q + P G = B''  sobre distintas mallas ===\n');
printf('  Se resuelve P y se mide cuánto se aparta de ser diagonal, si es\n');
printf('  positiva, y el residuo de la identidad con esa P.\n\n');
printf('  %-24s %-12s %-8s %-12s\n', 'malla', 'fuera diag', 'P>0', 'residuo');

casos = {};
for N = [10 40 160]
  casos{end+1} = {sprintf('uniforme N=%d', N), linspace(0, 1, N+1)};
end
for N = [10 40 160]
  t = linspace(0, 1, N+1);
  casos{end+1} = {sprintf('graduada t^2 N=%d', N), t.^2};
end
% La recuperación de P resuelve un sistema denso, de coste cúbico, así que
% las mallas adaptadas se comprueban hasta 800 celdas.
for N = [200 400 800]
  casos{end+1} = {sprintf('adaptada A N=%d', N), mallaAdaptada(N)};
end

for c = casos
  nom = c{1}{1};  x = c{1}{2};
  [G, D, B, ~] = operadoresMimeticos1D(x);
  Df = full(D);  Gf = full(G);  Bf = full(B);
  h  = diff(x(:).');
  Q  = diag([h(1), h, h(end)]);
  Popt = (Bf' - Df'*Q) / Gf;
  dp   = diag(Popt);
  offd = norm(Popt - diag(dp), inf);
  resi = norm(Df'*Q + diag(dp)*Gf - Bf', inf);
  printf('  %-24s %-12.2e %-8d %-12.2e\n', nom, offd, all(dp > 0), resi);
end

printf('\n=== Los pesos recuperados frente a Newton-Cotes ===\n');
printf('  Se muestra P_ii / hf_i. En el interior debe dar exactamente\n');
printf('  1, y 9/8 = 1.125 en la segunda cara y la penúltima.\n\n');
printf('  %-24s %-10s %-10s %-10s %-12s\n', 'malla', 'cara 1', 'cara 2', 'interior', 'cara N+1');

for c = casos
  nom = c{1}{1};  x = c{1}{2};
  [G, D, B, ~] = operadoresMimeticos1D(x);
  h  = diff(x(:).');
  Q  = diag([h(1), h, h(end)]);
  dp = diag((full(B)' - full(D)'*Q) / full(G));
  hf = [h(1), (h(1:end-1) + h(2:end))/2, h(end)];
  r  = dp(:).' ./ hf;
  printf('  %-24s %-10.5f %-10.5f %-10.5f %-12.5f\n', nom, r(1), r(2), ...
         max(abs(r(3:end-2) - 1)) + 1, r(end));
end
printf('\n');
end

% --------------------------------------------------------------------------
function x = mallaAdaptada(N)
  al = 250;  x0 = 0.8;
  k  = @(z) 1/al + al*(z-x0).^2;
  f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
  ga = al/(1+al^2*x0^2);
  gb = -atan(al*(1-x0)) - atan(al*x0);
  o  = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);
  x  = adaptaDeBoor(linspace(0, 1, N+1), k, f, ga, gb, o);
end
