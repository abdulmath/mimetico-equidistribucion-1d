function xnew = equidistribuye(x, rho, m, p, q)
% equidistribuye  Genera una malla que equidistribuye una función de
% densidad dada, iterando un número fijo de veces sobre h_i = C/rho_i.
%
%   xnew = equidistribuye(x, rho, m, p, q)
%
% Entrada
%   x   : caras de la malla actual, N+1 valores
%   rho : densidad de malla en los N centros de celda de x
%   m   : número de iteraciones internas (por omisión 3)
%   p,q : parámetros del suavizado de Dorfi y Drury
%
% En cada iteración la densidad se reevalúa sobre los centros de la malla en
% curso interpolando linealmente el perfil rho(x) de partida; sin ese paso la
% iteración sería idempotente y el lazo interno no tendría efecto alguno.
% Después se suaviza, se calcula C = (1/N) sum h_i rho_i, se toma h_i = C/rho_i
% y se renormaliza para que los h sumen la longitud del intervalo.

  if nargin < 3 || isempty(m), m = 3; end
  if nargin < 4, p = []; end
  if nargin < 5, q = []; end

  x   = x(:).';
  rho = rho(:).';
  N   = numel(x) - 1;
  a   = x(1);
  b   = x(end);

  xc0  = (x(1:N) + x(2:N+1))/2;      % centros de partida, soporte de rho
  xnew = x;

  for s = 1:m
    xc  = (xnew(1:N) + xnew(2:N+1))/2;
    r   = interp1(xc0, rho, xc, 'linear', 'extrap');
    r   = max(r, eps);
    r   = suavizaDorfiDrury(r, p, q);

    h   = diff(xnew);
    C   = sum(h .* r)/N;
    h   = C ./ r;
    h   = h * (b - a)/sum(h);

    xnew          = [a, a + cumsum(h)];
    xnew(N+1)     = b;
  end
end
