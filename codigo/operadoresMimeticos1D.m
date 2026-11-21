function [G, D, B, P, u, g, d] = operadoresMimeticos1D(x)
% operadoresMimeticos1D  Operadores miméticos de Castillo-Grone de orden 2
% sobre una malla escalonada no uniforme del intervalo [x(1), x(end)].
%
%   [G, D, B, P, u, g, d] = operadoresMimeticos1D(x)
%
% Entrada
%   x : vector fila de N+1 caras, estrictamente creciente.
%
% Salida
%   G : gradiente, de la malla mimética (N+2 valores) a las caras (N+1)
%   D : divergencia, de las caras (N+1) a la malla mimética (N+2)
%   B : matriz del teorema de Gauss discreto, (N+2) x (N+1)
%   P : pesos de Newton-Cotes sobre las caras, (N+1) x (N+1)
%   u : malla mimética [x(1), centros de celda, x(end)], N+2 valores
%   g : caras, donde vive el gradiente (coincide con x)
%   d : centros de celda, u(2:N+1)
%
% La primera y la última fila de G usan la plantilla de frontera (-8/3, 3,
% -1/3) escalada por la distancia s = -8/3 x_1 + 3 d_1 - 1/3 d_2, que sobre
% malla uniforme se reduce al paso h.

  x = x(:).';
  nnode = numel(x);
  N     = nnode - 1;

  if any(diff(x) <= 0)
    error('operadoresMimeticos1D: la malla debe ser estrictamente creciente');
  end

  u = [x(1), (x(1:N) + x(2:N+1))/2, x(end)];
  g = x;
  d = u(2:nnode);

  h = diff(x);

  % --- divergencia -------------------------------------------------------
  nucleo = spdiags([-ones(N,1), ones(N,1)], [0 1], N, N+1);
  D = sparse(N+2, N+1);
  D(2:N+1, :) = spdiags(1./h(:), 0, N, N) * nucleo;

  % --- gradiente ---------------------------------------------------------
  s      = zeros(1, nnode);
  s(1)   = -8/3*x(1) + 3*d(1) - 1/3*d(2);
  s(2:N) = d(2:N) - d(1:N-1);
  s(end) = 1/3*d(N-1) - 3*d(N) + 8/3*x(end);
  o      = 1./s;

  G = sparse(N+1, N+2);
  G(1, 1:3)         = o(1)   * [-8/3, 3, -1/3];
  G(N+1, N:N+2)     = o(end) * [ 1/3, -3, 8/3];
  for i = 2:N
    G(i, i)   = -o(i);
    G(i, i+1) =  o(i);
  end

  % --- matriz de Gauss ---------------------------------------------------
  B = sparse(N+2, N+1);
  B(1, 1)         = -1;
  B(2, 1:2)       = [ 1/8, -1/8];
  B(3, 1:2)       = [-1/8,  1/8];
  B(N,   N:N+1)   = [-1/8,  1/8];
  B(N+1, N:N+1)   = [ 1/8, -1/8];
  B(N+2, N+1)     =  1;

  % --- pesos de Newton-Cotes --------------------------------------------
  P = speye(N+1);
  P(1,1)     = 3/8;
  P(2,2)     = 9/8;
  P(N,N)     = 9/8;
  P(N+1,N+1) = 3/8;
end
