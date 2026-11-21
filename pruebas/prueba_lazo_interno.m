% prueba_lazo_interno.m
% Para una densidad fija, compara la recurrencia h_i = C/rho_i con la
% inversión exacta de la primitiva. La pregunta es si la recurrencia resuelve
% la misma ecuación de equidistribución cuando se la deja iterar, o si no
% converge a ella en absoluto.

function prueba_lazo_interno()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

a = 0;  b = 1;  alpha = 250;  x0 = 0.8;
gu = @(z) -(atan(alpha*(z-x0)) + atan(alpha*x0)) ...
          + (1-z).*(alpha./(1+alpha^2*(z-x0).^2));

N  = 100;
x  = linspace(a, b, N+1);
xc = (x(1:N) + x(2:N+1))/2;
r  = sqrt(1 + 0.5*gu(xc).^2);

xref = deBoor(x, r);

printf('\n=== Densidad fija, %d celdas ===\n', N);
printf('  m es el número de iteraciones internas de h = C/rho.\n');
printf('  %-5s %-15s %-15s\n', 'm', 'dispersión masa', 'dist. a de Boor');
for m = [1 2 3 5 10 20 50 100 200]
  xm = equidistribuye(x, r, m);
  printf('  %-5d %-15.3e %-15.3e\n', m, disper(masa(x, r, xm)), max(abs(xm - xref)));
end
printf('  %-5s %-15.3e %-15.3e\n', 'exacto', disper(masa(x, r, xref)), 0);
printf('\n');
end

% --------------------------------------------------------------------------
function d = disper(v)
  d = (max(v) - min(v))/mean(v);
end

function m = masa(x, rho, y)
  theta = [0, cumsum(diff(x).*rho)];
  t = interp1(x, theta, y, 'linear', 'extrap');
  m = diff(t);
end
