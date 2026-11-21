function [eInf, eL2] = erroresMalla(x, U, uex)
% erroresMalla  Error de la solución mimética en las normas del máximo y L2
% discreta sobre la malla escalonada.
%
%   [eInf, eL2] = erroresMalla(x, U, uex)
%
% Entrada
%   x   : caras de la malla, N+1 valores
%   U   : solución mimética, N+2 valores sobre [a, centros, b]
%   uex : manejador de la solución exacta
%
% La norma L2 se calcula como
%
%   ||e||_2^2 = sum_i h_i |e(x_{i+1/2})|^2 + (h_1/2)|e(a)|^2 + (h_N/2)|e(b)|^2,
%
% es decir, la regla del punto medio sobre cada celda más las dos medias
% celdas que aportan los valores de frontera. Los pesos son longitudes, de
% modo que la norma tiene las unidades correctas y no depende del número de
% nodos por un factor espurio.

  x = x(:).';
  U = U(:).';
  N = numel(x) - 1;
  u = [x(1), (x(1:N) + x(2:N+1))/2, x(end)];

  e = uex(u) - U;
  h = diff(x);

  eInf = max(abs(e));
  eL2  = sqrt(sum(h .* e(2:N+1).^2) + (h(1)/2)*e(1)^2 + (h(N)/2)*e(N+2)^2);
end
