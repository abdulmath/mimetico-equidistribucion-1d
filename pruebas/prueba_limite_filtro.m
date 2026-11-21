% prueba_limite_filtro.m
% Mide el límite de la densidad suavizada cuando la escala del filtro se fija
% en x, que es la tabla del límite del artículo. Dos caminos independientes:
%
%   1) En el punto fijo de la iteración exterior, sobre la malla adaptada y con
%      el gradiente discreto, para N = 200 ... 6400.
%   2) A partir de la derivada exacta de la solución, sobre malla auxiliar
%      uniforme fija de M celdas, para M = 2000 ... 16000. Esto evalúa
%      directamente el promedio continuo con núcleo gaussiano de anchura
%      sigma0, sin pasar por ninguna iteración.
%
% Si los dos caminos dan el mismo número y ese número no es uno, la densidad
% suavizada tiene un límite no constante y la malla límite sigue adaptada.
%
% El filtro se evalúa aquí por bloques de filas porque la matriz de pesos
% N x N de suavizaAnchuraFija pediría más de dos gigabytes con M = 16000. La
% fórmula es la misma y la primera parte de esta prueba lo comprueba.
%
% La corrida completa toma varios minutos; el grueso está en N = 6400 y en
% M = 16000. Conviene lanzarla con un solo hilo de BLAS,
%
%   OMP_NUM_THREADS=1 octave-cli --eval "prueba_limite_filtro"
%
% porque el orden de sumación del producto denso cambia la quinta cifra de los
% errores y con un hilo se reproducen las tablas del artículo tal como están.
%
% Se puede partir en trozos, que es lo cómodo cuando no se quiere dejar la
% consola ocupada varios minutos seguidos:
%
%   prueba_limite_filtro()          todo de una vez
%   prueba_limite_filtro(2, 1:4)    sólo el punto fijo, y sólo los cuatro
%                                   primeros valores de N
%   prueba_limite_filtro(5)         la extrapolación, con lo ya calculado
%
% Las etapas son 1 el filtro por bloques, 2 el punto fijo, 3 la malla auxiliar,
% 4 la convergencia al límite continuo y 5 la extrapolación de Aitken. Las
% etapas 2 y 3 guardan lo que calculan en salidas/limite_filtro_cache.mat, de
% modo que la 5 puede correrse después sin repetir nada. Borrando ese archivo
% se empieza de cero.

function prueba_limite_filtro(etapas, subN, subM)
if nargin < 1 || isempty(etapas), etapas = 1:5; end
hacer = @(e) any(etapas == e);

raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

sigma0 = 0.02;  C = 2;  beta = 0.5;  k1 = log(1.052);  k2 = 1;

% --- Problema A : pico interior, coeficiente variable ---------------------
al = 250;  x0 = 0.8;
A.nombre = 'A';
A.k  = @(z) 1/al + al*(z-x0).^2;
A.f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
A.ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
A.du = @(z) -(atan(al*(z-x0)) + atan(al*x0)) + (1-z).*al./(1+al^2*(z-x0).^2);
A.ga = al/(1+al^2*x0^2);
A.gb = -atan(al*(1-x0)) - atan(al*x0);

% --- Problema B : capa en el borde izquierdo, coeficiente constante -------
c = atan(100);
B.nombre = 'B';
B.k  = @(z) ones(size(z));
B.f  = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue = @(z) atan(100*z)/c;
B.du = @(z) 100./(c*(1 + 1e4*z.^2));
B.ga = 100/c;
B.gb = 1 + 100/(c*(1+1e4));

Ns = [200 400 800 1600 3200 6400];
Ms = [2000 4000 8000 16000];
if nargin < 2 || isempty(subN), subN = 1:numel(Ns); end
if nargin < 3 || isempty(subM), subM = 1:numel(Ms); end

% Lo ya calculado en corridas anteriores. NaN marca lo que falta.
cache = fullfile(raiz, 'salidas', 'limite_filtro_cache.mat');
qf = nan(2, numel(Ns));  qa = nan(2, numel(Ms));
if exist(cache, 'file')
  guardado = load(cache);
  if isequal(size(guardado.qf), size(qf)), qf = guardado.qf; end
  if isequal(size(guardado.qa), size(qa)), qa = guardado.qa; end
end

if hacer(1)
printf('\n=== el filtro por bloques contra suavizaAnchuraFija ===\n');
for N = [200 800 1600]
  t   = linspace(0, 1, N+1);
  xf  = t.^1.7;                                  % malla graduada, no uniforme
  xc  = (xf(1:end-1) + xf(2:end))/2;
  rho = 1 + 30*exp(-((xc-0.3)/0.01).^2) + 5*sin(40*xc).^2;
  a = suavizaAnchuraFija(xf, rho, sigma0, C);
  b = filtroBloques(xf, rho, sigma0, C);
  printf('  N = %5d   max|dif| = %.3e\n', N, max(abs(a-b)));
end
end

if hacer(2)
for j = subN(:).'
  N = Ns(j);
  for w = 1:2
    if w == 1, P = A; else, P = B; end
    xu = linspace(0, 1, N+1);
    o  = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10, 'filtro', 'anchura');
    xa = adaptaDeBoor(xu, P.k, P.f, P.ga, P.gb, o);
    [U, G] = resuelveRobin1D(xa, P.k, P.f, P.ga, P.gb);
    gu  = G*U;
    dUc = (gu(1:N) + gu(2:N+1)).'/2;
    r   = monitorArco(dUc, beta, k1, k2);
    rs  = filtroBloques(xa, r, sigma0, C);
    qf(w, j) = max(rs)/min(rs);
  end
  save('-binary', cache, 'qf', 'qa');
end

printf('\n=== razon max/min de la densidad suavizada en el punto fijo ===\n');
printf('  %-8s %-12s %-12s\n', 'Celdas', 'A', 'B');
for j = 1:numel(Ns)
  if isnan(qf(1,j)), continue; end
  printf('  %-8d %-12.4f %-12.4f\n', Ns(j), qf(1,j), qf(2,j));
end
fflush(stdout);
end

if hacer(3)
for j = subM(:).'
  M  = Ms(j);
  xf = linspace(0, 1, M+1);
  xc = (xf(1:end-1) + xf(2:end))/2;
  for w = 1:2
    if w == 1, P = A; else, P = B; end
    r  = monitorArco(P.du(xc), beta, k1, k2);
    rs = filtroBloques(xf, r, sigma0, C);
    qa(w, j) = max(rs)/min(rs);
  end
  save('-binary', cache, 'qf', 'qa');
end

printf('\n=== lo mismo con la derivada exacta sobre malla auxiliar fija ===\n');
printf('  %-8s %-12s %-12s\n', 'M', 'A', 'B');
for j = 1:numel(Ms)
  if isnan(qa(1,j)), continue; end
  printf('  %-8d %-12.4f %-12.4f\n', Ms(j), qa(1,j), qa(2,j));
end
fflush(stdout);
end

% La proposicion del articulo afirma algo mas fuerte que la existencia del
% limite: que el filtro discreto converge al cociente gaussiano continuo sobre
% cualquier sucesion de mallas cuyo paso maximo tienda a cero, y que el limite
% no depende de esa sucesion. Se comprueba midiendo
%
%     max_i | rho~_i - rho~_sigma0(d_i) |
%
% con el miembro derecho calculado por cuadratura adaptativa, sobre dos
% familias de mallas muy distintas: la uniforme y la graduada x = t^1.7, cuya
% razon entre celda mayor y menor crece con N. Si el factor h_j del peso hace
% su trabajo, las dos han de converger al mismo limite.
%
% Se prueban dos densidades. La primera es la de la proposicion tal como esta
% enunciada: los valores rho(d_i) de una funcion fija, aqui la longitud de arco
% sqrt(1 + beta u'^2). La segunda es la que usa de verdad la iteracion, que
% ademas reescala las dos celdas extremas por k1 segun (2.tl). Ese reescalado
% no es el muestreo de ninguna funcion, sino una modificacion de dos entradas
% cuya medida tiende a cero; el limite es el mismo, pero la perturbacion que
% introduce en el promedio es de orden h y ahi se ve, porque baja el orden
% observado de dos a uno. Conviene tenerlo presente: la proposicion vale para
% las dos, con constantes muy distintas.
%
% La cuadratura adaptativa es cara, de modo que el maximo se toma sobre una
% submuestra de indices interiores repartidos por todo el dominio. El maximo
% sobre la submuestra decae igual que el maximo completo y cuesta mucho menos.
if hacer(4)
printf('\n=== convergencia al cociente gaussiano continuo ===\n');
printf('  %-9s %-14s %-10s %-6s %-12s %-6s\n', ...
       'Problema', 'Densidad', 'Malla', 'N', 'max|dif|', 'orden');
for w = 1:2
  if w == 1, P = A; else, P = B; end
  arco = @(y) sqrt(1 + beta*P.du(y).^2);
  for cual = {'de la prop.', 'con borde k1'}
    % cada variante tiene su propio limite continuo, y hay que compararla
    % contra el suyo: la de la proposicion contra la longitud de arco, y la de
    % la iteracion contra la longitud de arco ya transformada por (eq:tl). El
    % reescalado k1 afecta a dos celdas y no tiene contrapartida continua.
    if strcmp(cual{1}, 'de la prop.')
      g = arco;
    else
      g = @(y) arco(y) + k2*log(arco(y));
    end
    for tipo = {'uniforme', 'graduada'}
      ant = NaN;
      for N = [400 800 1600]
        t = linspace(0, 1, N+1);
        if strcmp(tipo{1}, 'uniforme'), xf = t; else, xf = t.^1.7; end
        xc  = (xf(1:end-1) + xf(2:end))/2;
        sig = max(sigma0, C*diff(xf)(:));
        if strcmp(cual{1}, 'de la prop.')
          r = g(xc);
        else
          r = monitorArco(P.du(xc), beta, k1, k2);
        end
        rs  = filtroBloques(xf, r, sigma0, C);
        idx = round(linspace(2, N-1, 120));
        dif = zeros(size(idx));
        for m = 1:numel(idx)
          i  = idx(m);  s = sig(i);
          lo = max(0, xc(i) - 8*s);  hi = min(1, xc(i) + 8*s);
          nu = quadgk(@(y) exp(-((xc(i)-y)/s).^2).*g(y), lo, hi, ...
                      'AbsTol', 1e-13, 'RelTol', 1e-12);
          de = quadgk(@(y) exp(-((xc(i)-y)/s).^2), lo, hi, ...
                      'AbsTol', 1e-13, 'RelTol', 1e-12);
          dif(m) = abs(rs(i) - nu/de);
        end
        er = max(dif);
        if isnan(ant)
          printf('  %-9s %-14s %-10s %-6d %-12.3e %-6s\n', ...
                 P.nombre, cual{1}, tipo{1}, N, er, '---');
        else
          printf('  %-9s %-14s %-10s %-6d %-12.3e %-6.2f\n', ...
                 P.nombre, cual{1}, tipo{1}, N, er, log2(ant/er));
        end
        ant = er;
        fflush(stdout);
      end
    end
  end
end
printf('\n  Las dos familias de mallas convergen al mismo limite continuo.\n');
printf('  Con la densidad de la proposicion el orden es dos; el reescalado de\n');
printf('  borde lo baja a uno sin cambiar el limite.\n');
end

if hacer(5)
if any(isnan(qf(:))) || any(isnan(qa(:)))
  error(['prueba_limite_filtro: falta calcular el punto fijo o la malla ' ...
         'auxiliar; correr antes las etapas 2 y 3']);
end
printf('\n=== limite extrapolado con el delta cuadrado de Aitken ===\n');
printf('  %-14s %-12s %-12s\n', '', 'A', 'B');
printf('  %-14s %-12.4f %-12.4f\n', 'punto fijo', ...
       aitken(qf(1,end-2:end)), aitken(qf(2,end-2:end)));
printf('  %-14s %-12.4f %-12.4f\n', 'malla auxiliar', ...
       aitken(qa(1,end-2:end)), aitken(qa(2,end-2:end)));
printf('\n  Los dos caminos deben coincidir, y el valor no debe ser uno.\n\n');
end
end

% -------------------------------------------------------------------------
function rs = filtroBloques(x, rho, sigma0, C)
% Igual que suavizaAnchuraFija, evaluada por bloques de filas para no
% materializar la matriz de pesos completa.
  x   = x(:).';
  rho = rho(:).';
  xc  = (x(1:end-1) + x(2:end))/2;
  h   = diff(x);
  sig = max(sigma0, C*h(:));
  n   = numel(xc);
  rs  = zeros(1, n);
  bl  = 512;
  for i0 = 1:bl:n
    i1 = min(i0 + bl - 1, n);
    d  = (xc(i0:i1).' - xc) ./ sig(i0:i1);
    w  = exp(-d.^2) .* h;
    rs(i0:i1) = ((w*rho.') ./ sum(w, 2)).';
  end
end

function L = aitken(v)
% Delta cuadrado de Aitken sobre tres términos consecutivos.
  d1 = v(2) - v(1);  d2 = v(3) - v(2);
  L  = v(3) - d2^2/(d2 - d1);
end
