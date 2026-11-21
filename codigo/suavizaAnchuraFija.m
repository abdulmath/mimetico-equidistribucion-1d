function rs = suavizaAnchuraFija(x, rho, sigma0, C)
% suavizaAnchuraFija  Filtro de paso bajo sobre la densidad de malla cuya
% ventana tiene anchura fija en la variable x, no un número fijo de celdas.
%
%   rs = suavizaAnchuraFija(x, rho, sigma0, C)
%
% Entrada
%   x      : caras de la malla, N+1 valores crecientes
%   rho    : densidad en los centros de celda, N valores
%   sigma0 : anchura mínima de la ventana, en unidades de x
%   C      : número de celdas que fija la anchura mientras la malla es gruesa
%
% Cada valor se sustituye por el promedio de rho con núcleo gaussiano
%
%   w_ij = exp(-((xc_i - xc_j)/sigma_i)^2) h_j,     sigma_i = max(sigma0, C h_i)
%
% El peso h_j hace que el promedio no dependa de cómo estén repartidas las
% celdas, que es justo lo que se busca: sobre una malla muy poco uniforme el
% resultado debe ser el mismo que sobre una uniforme.
%
% El motivo de esta forma de suavizar es que el filtro clásico de Dorfi y
% Drury promedia un número fijo de celdas, de modo que su ventana se encoge
% al refinar y la densidad nunca llega a ser suave como función de x. La malla
% que resulta conserva estructura a escala de celda, la razón h_{i+1}/h_i no
% tiende a uno, y el esquema mimético pierde la supraconvergencia. Con anchura
% fija en x esa razón decae como O(1/N) y el orden dos se recupera.
%
% El término C h_i sólo actúa mientras h_i > sigma0/C. Sirve para que en las
% mallas gruesas la ventana abarque al menos unas cuantas celdas y el punto
% fijo exterior siga siendo atractivo; en cuanto la malla se refina, sigma0
% domina y el filtro pasa a tener anchura constante.
%
% Valores por omisión: sigma0 = 0.02 sobre [0,1] reescalado, C = 2.

  if nargin < 3 || isempty(sigma0), sigma0 = 0.02*(x(end) - x(1)); end
  if nargin < 4 || isempty(C),      C      = 2;                    end

  x   = x(:).';
  rho = rho(:).';
  xc  = (x(1:end-1) + x(2:end))/2;
  h   = diff(x);

  sig = max(sigma0, C*h(:));
  d   = (xc.' - xc) ./ sig;
  w   = exp(-d.^2) .* h;
  rs  = ((w*rho.') ./ sum(w, 2)).';
end
