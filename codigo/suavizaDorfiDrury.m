function rs = suavizaDorfiDrury(rho, p, q)
% suavizaDorfiDrury  Filtro de paso bajo de Dorfi y Drury sobre la función
% de densidad de malla.
%
%   rs = suavizaDorfiDrury(rho, p, q)
%
% Cada valor se sustituye por el promedio ponderado de sus 2p+1 vecinos con
% pesos (q/(q+1))^{|k-i|}. Cerca de los extremos la ventana se recorta y los
% pesos se renormalizan, de modo que rs conserva el rango de rho.
%
% Valores por omisión: p = 3, q = 2.

  if nargin < 2 || isempty(p), p = 3; end
  if nargin < 3 || isempty(q), q = 2; end

  rho = rho(:).';
  n   = numel(rho);
  w   = (q/(q+1)).^(0:p);
  rs  = zeros(1, n);

  for i = 1:n
    lo = max(1, i-p);
    hi = min(n, i+p);
    ww = w(abs((lo:hi) - i) + 1);
    rs(i) = (ww * rho(lo:hi).') / sum(ww);
  end
end
