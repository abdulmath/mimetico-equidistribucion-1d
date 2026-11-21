function rho = monitorArco(dU, beta, k1, k2)
% monitorArco  Función de densidad de malla de longitud de arco escalada,
% opcionalmente transformada según Tapia y López.
%
%   rho = monitorArco(dU, beta)              longitud de arco escalada
%   rho = monitorArco(dU, beta, k1, k2)      además la transformación
%
% Entrada
%   dU   : derivada aproximada de la solución en los centros de celda
%   beta : parámetro de escala, positivo
%   k1   : factor de los dos valores extremos
%   k2   : peso del término logarítmico en los valores interiores
%
% La transformación refuerza la densidad donde la curvatura es grande; se
% aplica antes del suavizado, no después.

  dU  = dU(:).';
  rho = sqrt(1 + beta*dU.^2);

  if nargin >= 3 && ~isempty(k1)
    n = numel(rho);
    r = rho;
    r(1)       = k1*rho(1);
    r(n)       = k1*rho(n);
    r(2:n-1)   = rho(2:n-1) + k2*log(rho(2:n-1));
    rho        = max(r, eps);
  end
end
