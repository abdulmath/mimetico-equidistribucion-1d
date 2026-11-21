function [x, U, hist] = adaptaDeBoor(x, k, f, ga, gb, opts)
% adaptaDeBoor  Malla adaptada en movimiento donde cada paso equidistribuye
% exactamente la densidad en curso mediante inversión de su primitiva.
%
%   [x, U, hist] = adaptaDeBoor(x, k, f, ga, gb, opts)
%
% El lazo exterior es un punto fijo sobre la malla: se resuelve, se calcula
% la densidad a partir del gradiente mimético, se suaviza y se equidistribuye.
% Lo que hace contractiva la iteración es que el paso de equidistribución sea
% exacto. El promedio
%
%   x^{j+1} = (1-omega) x^j + omega deBoor(x^j, rho^j)
%
% no cambia el punto fijo, sólo la rapidez con que se llega a él. Con omega=1
% el lazo converge en menos de quince pasos siempre que la malla resuelva la
% capa; por debajo de unas doscientas celdas para una capa de anchura 1/250 el
% lazo entra en un ciclo y hace falta omega en torno a 0.7, que converge en
% los dos problemas de prueba y a la misma malla, con cinco cifras iguales.
%
% La densidad se suaviza con un filtro cuya ventana tiene anchura fija en x.
% Esto no es un detalle de ajuste. El filtro clásico de Dorfi y Drury promedia
% un número fijo de celdas, así que su ventana se encoge al refinar, la
% densidad conserva estructura a escala de celda y la malla resultante no
% alcanza la suavidad h_{i+1}/h_i = 1 + O(h) que pide la supraconvergencia del
% esquema mimético. Con anchura fija en x la razón de tamaños decae como
% O(1/N) y el orden dos se recupera. El filtro de Dorfi y Drury sigue
% disponible mediante opts.filtro para poder comparar.
%
% Campos de opts (todos opcionales)
%   beta   parámetro de escala del monitor                (0.5)
%   k1,k2  parámetros de la transformación de Tapia-López (log(1.052), 1)
%   omega  factor de relajación en (0,1]                  (1)
%   filtro 'anchura', 'dorfidrury' o 'global'             ('anchura')
%   betaS  parametro del suavizado global en xi           (20)
%   sigma0 anchura mínima del filtro, en unidades de x    (0.02*(b-a))
%   C      celdas que fijan la anchura en malla gruesa    (2)
%   p,q    parámetros del suavizado de Dorfi y Drury      (3, 2)
%   maxIt  tope de iteraciones exteriores                 (200)
%   tolX   tolerancia sobre el movimiento de la malla     (1e-10)

  if nargin < 6, opts = struct(); end
  beta   = valor(opts, 'beta',   0.5);
  k1     = valor(opts, 'k1',     log(1.052));
  k2     = valor(opts, 'k2',     1);
  omega  = valor(opts, 'omega',  1);
  filtro = valor(opts, 'filtro', 'anchura');
  sigma0 = valor(opts, 'sigma0', 0.02*(x(end) - x(1)));
  C      = valor(opts, 'C',      2);
  p      = valor(opts, 'p',      3);
  q      = valor(opts, 'q',      2);
  betaS  = valor(opts, 'betaS',  20);
  maxIt  = valor(opts, 'maxIt',  200);
  tolX   = valor(opts, 'tolX',   1e-10);

  x = x(:).';
  N = numel(x) - 1;
  hist = struct('mov', [], 'it', 0, 'convergio', false);

  for j = 1:maxIt
    [U, G] = resuelveRobin1D(x, k, f, ga, gb);

    gu  = G*U;
    dUc = (gu(1:N) + gu(2:N+1)).'/2;

    rho = monitorArco(dUc, beta, k1, k2);
    switch filtro
      case 'dorfidrury'
        rho = suavizaDorfiDrury(rho, p, q);
      case 'global'
        rho = suavizaGlobal(rho, betaS);
      otherwise
        rho = suavizaAnchuraFija(x, rho, sigma0, C);
    end

    xeq = deBoor(x, rho);
    xn  = (1-omega)*x + omega*xeq;
    xn(1) = x(1);  xn(N+1) = x(end);

    mov = max(abs(xn - x));
    hist.mov(end+1) = mov;
    hist.it = j;
    x = xn;

    if mov < tolX
      hist.convergio = true;
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
