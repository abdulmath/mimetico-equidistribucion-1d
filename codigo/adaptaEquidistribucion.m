function [x, U, hist] = adaptaEquidistribucion(x, k, f, ga, gb, opts)
% adaptaEquidistribucion  Malla adaptada en movimiento por equidistribución
% de la longitud de arco escalada, con continuación sobre los parámetros del
% monitor.
%
%   [x, U, hist] = adaptaEquidistribucion(x, k, f, ga, gb, opts)
%
% El número de nodos permanece fijo: sólo se mueven. En cada iteración
% exterior se resuelve el problema en la malla en curso, se construye la
% densidad de malla a partir del gradiente mimético, se transforma según
% Tapia y López, y se generan las nuevas posiciones con equidistribuye.
%
% La continuación sobre el parámetro de escala beta y sobre el peso k2 es la
% que evita que la iteración exterior entre en un ciclo: se arranca con un
% monitor casi plano y se refuerza gradualmente, de modo que cada malla es
% una perturbación pequeña de la anterior.
%
% Campos de opts (todos opcionales)
%   beta   valor final del parámetro de escala           (0.5)
%   k1     factor de los extremos en Tapia-López         (log(1.052))
%   k2     valor final del peso logarítmico              (1)
%   nCont  iteraciones de la rampa de continuación       (10)
%   mInt   iteraciones internas de equidistribuye        (3)
%   p, q   parámetros del suavizado                      (3, 2)
%   maxIt  iteraciones exteriores como tope              (60)
%   tolX   tolerancia sobre el movimiento de la malla    (1e-8)
%
% hist devuelve, por iteración, el movimiento máximo de la malla, para poder
% comprobar la convergencia del proceso.

  if nargin < 6, opts = struct(); end
  beta  = valor(opts, 'beta',  0.5);
  k1    = valor(opts, 'k1',    log(1.052));
  k2    = valor(opts, 'k2',    1);
  nCont = valor(opts, 'nCont', 10);
  mInt  = valor(opts, 'mInt',  3);
  p     = valor(opts, 'p',     3);
  q     = valor(opts, 'q',     2);
  maxIt = valor(opts, 'maxIt', 60);
  tolX  = valor(opts, 'tolX',  1e-8);

  x = x(:).';
  N = numel(x) - 1;
  hist = struct('mov', [], 'it', 0);

  for j = 1:maxIt
    [U, G] = resuelveRobin1D(x, k, f, ga, gb);

    gu  = G*U;                                   % gradiente en las caras
    dUc = (gu(1:N) + gu(2:N+1)).'/2;             % llevado a los centros

    t  = min(1, j/nCont);                        % rampa de continuación
    bj = beta*t;
    kj = k2*t;

    rho  = monitorArco(dUc, bj, k1, kj);
    xnew = equidistribuye(x, rho, mInt, p, q);

    mov = max(abs(xnew - x));
    hist.mov(end+1) = mov;
    hist.it = j;
    x = xnew;

    if j > nCont && mov < tolX
      break
    end
  end

  U = resuelveRobin1D(x, k, f, ga, gb);
end

% --------------------------------------------------------------------------
function v = valor(s, campo, porOmision)
  if isfield(s, campo) && ~isempty(s.(campo))
    v = s.(campo);
  else
    v = porOmision;
  end
end
