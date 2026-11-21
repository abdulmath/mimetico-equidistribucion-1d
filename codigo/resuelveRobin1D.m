function [U, G, u, g, d] = resuelveRobin1D(x, k, f, ga, gb)
% resuelveRobin1D  Resuelve (k(x) u')' = f(x) en [x(1), x(end)] con
% condiciones de Robin u(a) + u'(a) = ga, u(b) + u'(b) = gb, mediante el
% esquema mimético de Castillo-Grone de orden 2 sobre malla no uniforme.
%
%   [U, G, u, g, d] = resuelveRobin1D(x, k, f, ga, gb)
%
% Entrada
%   x  : caras de la malla, vector fila creciente
%   k  : manejador de función, coeficiente evaluado en las caras
%   f  : manejador de función, término fuente evaluado en los centros
%   ga : valor de u(a) + u'(a)
%   gb : valor de u(b) + u'(b)
%
% Salida
%   U : solución sobre la malla mimética, vector columna de N+2 valores
%   G : operador gradiente empleado, por si hace falta el flujo o el monitor

  [G, D, B, P, u, g, d] = operadoresMimeticos1D(x);
  N = numel(x) - 1;

  A = sparse([1, N+2], [1, N+2], [-1, 1], N+2, N+2);
  K = spdiags(k(g(:)), 0, N+1, N+1);

  M = A + B*G + D*K*G;
  F = [-ga; f(d(:)); gb];

  U = M \ F;
end
