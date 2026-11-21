function rs = suavizaGlobal(rho, beta)
% suavizaGlobal  Suavizado global de la funcion monitor resolviendo
%
%   (I - beta^{-2} d^2/dxi^2) rs = rho
%
% en la variable computacional xi, con condiciones de Neumann discretas.
% Es el esquema (2.132) de Huang y Russell, del que el filtro local de pesos
% geometricos es la truncacion. Con beta fijo, la cota del teorema 2.5.2
%
%   nu <= rs_j/rs_{j-1} <= 1/nu,   nu = (sqrt(1+4/(beta^2 dxi^2)) - 1)
%                                     /(sqrt(1+4/(beta^2 dxi^2)) + 1)
%
% tiende a uno al refinar, porque dxi = 1/N -> 0.

  rho = rho(:);
  n   = numel(rho);
  dxi = 1/n;
  a   = 1/(beta^2*dxi^2);

  sub  = [-a*ones(n-2,1); -1; 0];        % subdiagonal, indices 1..n-1 usados
  sup  = [0; -1; -a*ones(n-2,1)];        % superdiagonal
  dia  = [1; (1+2*a)*ones(n-2,1); 1];
  A    = spdiags([sub, dia, sup], -1:1, n, n);

  b       = rho;
  b(1)    = 0;
  b(n)    = 0;

  rs = (A\b).';
end
