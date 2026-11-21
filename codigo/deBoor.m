function xnew = deBoor(x, rho, M)
% deBoor  Malla que equidistribuye exactamente una densidad constante a
% trozos, por inversión de su primitiva.
%
%   xnew = deBoor(x, rho)      malla nueva con el mismo número de celdas
%   xnew = deBoor(x, rho, M)   malla nueva con M celdas
%
% Entrada
%   x   : caras de la malla sobre la que está dada la densidad, N+1 valores
%   rho : densidad, un valor constante por celda, N valores
%
% Con rho constante en cada celda la primitiva
%
%   theta(z) = int_a^z rho(s) ds
%
% es continua y lineal a trozos, luego estrictamente creciente e invertible
% de forma exacta. La malla equidistribuida es
%
%   xnew_j = theta^{-1}( (j-1) theta(b) / M ),   j = 1, ..., M+1,
%
% y no hace falta iterar: el resultado equidistribuye la densidad dada hasta
% el error de redondeo. Iterar sólo tiene sentido cuando la propia densidad
% depende de la malla, que es lo que hace el lazo exterior.

  x   = x(:).';
  rho = rho(:).';
  N   = numel(x) - 1;
  if nargin < 3 || isempty(M), M = N; end

  if numel(rho) ~= N
    error('deBoor: rho debe tener un valor por celda');
  end
  if any(rho <= 0)
    error('deBoor: la densidad debe ser positiva');
  end

  h     = diff(x);
  theta = [0, cumsum(h .* rho)];      % primitiva en las caras
  total = theta(end);

  objetivo = (0:M) * total / M;

  xnew    = zeros(1, M+1);
  xnew(1) = x(1);
  xnew(M+1) = x(end);

  i = 1;
  for j = 2:M
    t = objetivo(j);
    while i < N && theta(i+1) < t
      i = i + 1;
    end
    xnew(j) = x(i) + (t - theta(i))/rho(i);
  end
end
