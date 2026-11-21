% figuras_articulo.m
% Genera las figuras del artículo. Cada una responde a una afirmación concreta
% del texto y todas se recalculan aquí: ningún número está copiado a mano.
%
% Cada panel se escribe en su propio archivo, para poder editarlo por separado
% y componerlo después en LaTeX con subfiguras. Los nombres son
%
%   fig_solucion_A      solución exacta y solución mimética sobre la adaptada
%   fig_solucion_B
%   fig_caras_A         caras de la malla uniforme y de la adaptada
%   fig_caras_B
%   fig_tamano_A        tamaño de celda frente a x, en escala logarítmica
%   fig_tamano_B
%   fig_salto_A         salto máximo |h_{i+1}/h_i - 1| frente a N con los dos
%   fig_salto_B         filtros; es la figura del resultado central
%   fig_convergencia_A  error frente a número de celdas, uniforme y adaptada
%   fig_convergencia_B
%
% De cada uno salen dos archivos, .pdf y .png. El PDF es vectorial: es el que
% entra en el artículo, se amplía sin perder nitidez y se puede abrir en un
% editor de dibujo para retocar rótulos, trazos y colores. El PNG queda para
% previsualizar y para donde no se acepte vectorial. Añadiendo '-dsvg' a
% FORMATOS sale además el SVG, que es el formato cómodo para Inkscape.
%
% Los paneles no llevan título dentro del dibujo: el nombre del problema va en
% la subleyenda de LaTeX, que es donde le corresponde.
%
% Las definiciones de los dos problemas son las mismas de prueba_dos_problemas,
% de manera que las figuras y las tablas del artículo salen del mismo cálculo.
%
% Desde la consola:  octave-cli --eval "figuras_articulo"
%
% Admite un argumento opcional con las etapas que se quieren rehacer, por si
% sólo hace falta retocar una figura:
%
%   figuras_articulo()        las cuatro etapas
%   figuras_articulo([3 4])   sólo el salto y la convergencia
%
% Las etapas son 1 soluciones y caras, 2 tamaño de celda, 3 salto máximo,
% 4 convergencia.
%
% Un segundo argumento elige el idioma de los rótulos de dentro del dibujo:
%
%   figuras_articulo([], 'es')   rótulos en español, archivos sin sufijo
%   figuras_articulo([], 'en')   rótulos en inglés, archivos con sufijo _en
%
% Las dos versiones conviven en figuras/ y cada artículo incluye la suya. Los
% rótulos que son los mismos en los dos idiomas, como x, N, u(x) y h_i, no se
% traducen.

function figuras_articulo(etapas, idioma)

if nargin < 1 || isempty(etapas), etapas = 1:4; end
if nargin < 2 || isempty(idioma), idioma = 'es'; end
hacer = @(e) any(etapas == e);

% Rótulos de dentro del dibujo. Con idioma = 'en' los paneles se escriben con
% el sufijo _en, de modo que las dos versiones conviven en figuras/ y cada
% artículo incluye la suya.
switch lower(idioma)
  case 'en'
    R = struct('exacta','exact', 'adaptada','adapted', 'uniforme','uniform', ...
               'anchura','fixed width', 'dd','Dorfi and Drury', ...
               'pend1','slope -1', 'pend2','slope -2', ...
               'err','maximum norm error');
    SUF = '_en';
  case 'es'
    R = struct('exacta','exacta', 'adaptada','adaptada', 'uniforme','uniforme', ...
               'anchura','anchura fija', 'dd','Dorfi y Drury', ...
               'pend1','pendiente -1', 'pend2','pendiente -2', ...
               'err','error en norma del máximo');
    SUF = '';
  otherwise
    error('figuras_articulo: idioma debe ser ''es'' o ''en''');
end

global alpha x0

raiz   = fileparts(fileparts(mfilename('fullpath')));
dirFig = fullfile(raiz, 'figuras');
addpath(fullfile(raiz, 'codigo'));
if ~exist(dirFig, 'dir'), mkdir(dirFig); end
backendGrafico();

FORMATOS = {'-dpdf', '-dpng'};   % añádase '-dsvg' si se va a editar en Inkscape

% El panel se dibuja al tamaño con el que va a imprimirse, no más grande. La
% caja de texto del artículo mide 15.8 cm, o sea 449 puntos, y cada subfigura
% ocupa 0.49 de ella: 220 puntos. Dibujando a ese ancho, LaTeX no reescala y
% los rótulos salen del tamaño que se les pone aquí. Si se dibujara al doble y
% se redujera después, la letra saldría a la mitad y no habría quien la lea.
ANCHO    = 220;
ALTO     = 165;
LETRA    = 8;                    % puntos, contra los 11 del cuerpo del texto

% --- los dos problemas de prueba -----------------------------------------
al = 250;  xp = 0.8;  c = atan(100);
alpha = al;  x0 = xp;

A.nom = 'Problema A: pico interior en x_0 = 0.8';
A.let = 'A';
A.k   = @(z) 1/al + al*(z-xp).^2;
A.f   = @(z) -2*(1 + al*(z-xp).*(atan(al*xp) + atan(al*(z-xp))));
A.ue  = @(z) (1-z).*(atan(al*(z-xp)) + atan(al*xp));
A.gu  = @(z) -(atan(al*(z-xp)) + atan(al*xp)) + (1-z).*(al./(1+al^2*(z-xp).^2));
A.ga  = al/(1+al^2*xp^2);
A.gb  = -atan(al*(1-xp)) - atan(al*xp);
A.cf  = [-al/(1+al^2*xp^2), -atan(al*(1-xp))-atan(al*xp), -1, 1];
A.leg = 'northwest';   % el pico está a la derecha, la esquina libre es esa

B.nom = 'Problema B: capa en x = 0';
B.let = 'B';
B.k   = @(z) ones(size(z));
B.f   = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue  = @(z) atan(100*z)/c;
B.gu  = @(z) 100./(c*(1 + 1e4*z.^2));
B.ga  = 100/c;
B.gb  = 1 + 100/(c*(1+1e4));
B.cf  = [-100/c, 1 + 100/((1+1e4)*c), 1, 1];
B.leg = 'southeast';   % la capa sube en x = 0 y la curva se queda arriba

P  = [A B];
oA = struct('omega', 1, 'maxIt', 300, 'tolX', 1e-10);
oD = oA;  oD.filtro = 'dorfidrury';
oM = oA;  oM.omega  = 0.7;
Nf = 200;    % malla de las figuras del tamaño de celda
Nm = 100;    % malla del panel de caras, donde importa poder contarlas. Por
             % debajo de unas cien celdas la ventana del filtro se traga la
             % capa del problema B y la malla sale casi uniforme.
Ns = [200 400 800 1600 3200];

% =========================================================================
% Las soluciones y las caras de las mallas
% =========================================================================
if hacer(1)
printf('=== 1. Soluciones y mallas adaptadas ===\n');  fflush(stdout);

for s = 1:2
  p  = P(s);
  xu = linspace(0, 1, Nf+1);
  [xa, Ua] = adaptaDeBoor(xu, p.k, p.f, p.ga, p.gb, oA);
  ua = [xa(1), (xa(1:end-1)+xa(2:end))/2, xa(end)];
  zz = linspace(0, 1, 4000);

  f = nuevaFigura(ANCHO, ALTO, LETRA);
  plot(zz, p.ue(zz), 'k-', 'linewidth', 1.4); hold on;
  plot(ua, Ua, 'ro', 'markersize', 2.5);
  xlabel('x'); ylabel('u(x)'); grid on;
  leyenda(LETRA-1, p.leg, R.exacta, R.adaptada);
  guarda(f, dirFig, ['fig_solucion_' p.let SUF], FORMATOS);

  % El panel de caras usa una malla gruesa para que se puedan contar. Por
  % debajo de unas doscientas celdas el lazo con omega = 1 entra en ciclo,
  % así que aquí se relaja, que es lo que dice la sección del lazo exterior.
  [xm, ~, hm] = adaptaDeBoor(linspace(0, 1, Nm+1), p.k, p.f, p.ga, p.gb, oM);
  xv = linspace(0, 1, Nm+1);
  printf('  malla del panel: %d celdas, %d iteraciones, converge %d\n', ...
         Nm, hm.it, hm.convergio);

  f = nuevaFigura(ANCHO, round(ALTO*0.62), LETRA);
  plot([xm; xm], repmat([0; 0.42], 1, numel(xm)), 'b-', 'linewidth', 0.7);
  hold on;
  plot([xv; xv], repmat([1; 1.42], 1, numel(xv)), 'k-', 'linewidth', 0.7);
  set(gca, 'ytick', [0.21 1.21], 'yticklabel', {R.adaptada, R.uniforme});
  xlim([0 1]); ylim([-0.25 1.7]); xlabel('x');
  guarda(f, dirFig, ['fig_caras_' p.let SUF], FORMATOS);

  printf('  %s: h_max/h_min = %.1f\n', p.nom, max(diff(xa))/min(diff(xa)));
  fflush(stdout);
end
end

% =========================================================================
% El tamaño de celda frente a x
% =========================================================================
if hacer(2)
printf('=== 2. Tamaño de celda ===\n');  fflush(stdout);

for s = 1:2
  p  = P(s);
  xu = linspace(0, 1, Nf+1);
  xa = adaptaDeBoor(xu, p.k, p.f, p.ga, p.gb, oA);
  xd = adaptaDeBoor(xu, p.k, p.f, p.ga, p.gb, oD);

  f = nuevaFigura(ANCHO, ALTO, LETRA);
  semilogy(centros(xa), diff(xa), 'b-', 'linewidth', 1.3); hold on;
  semilogy(centros(xd), diff(xd), 'r--', 'linewidth', 1.1);
  semilogy([0 1], (1/Nf)*[1 1], 'k:', 'linewidth', 1.2);
  xlabel('x'); ylabel('h_i'); grid on;
  leyenda(LETRA-1, 'southeast', R.anchura, R.dd, R.uniforme);
  guarda(f, dirFig, ['fig_tamano_' p.let SUF], FORMATOS);
end
end

% =========================================================================
% El salto máximo frente a N
% =========================================================================
if hacer(3)
printf('=== 3. Salto máximo de la razón de tamaños ===\n');  fflush(stdout);

Sa = zeros(numel(Ns), 2);  Sd = zeros(numel(Ns), 2);
for s = 1:2
  p = P(s);
  for m = 1:numel(Ns)
    xu = linspace(0, 1, Ns(m)+1);
    xa = adaptaDeBoor(xu, p.k, p.f, p.ga, p.gb, oA);
    xd = adaptaDeBoor(xu, p.k, p.f, p.ga, p.gb, oD);
    Sa(m,s) = salto(xa);  Sd(m,s) = salto(xd);
    printf('  %s  N=%4d   anchura %.4f   Dorfi-Drury %.4f\n', ...
           p.let, Ns(m), Sa(m,s), Sd(m,s));
    fflush(stdout);
  end
end

for s = 1:2
  f = nuevaFigura(ANCHO, ALTO, LETRA);
  loglog(Ns, Sa(:,s), 'bo-', 'linewidth', 1.4, 'markersize', 5); hold on;
  loglog(Ns, Sd(:,s), 'rs--', 'linewidth', 1.2, 'markersize', 5);
  loglog(Ns, Sa(1,s)*Ns(1)./Ns, 'k:', 'linewidth', 1.2);
  xlabel('N'); ylabel('max_i |h_{i+1}/h_i - 1|'); grid on;
  leyenda(LETRA-1, 'southwest', R.anchura, R.dd, R.pend1);
  guarda(f, dirFig, ['fig_salto_' P(s).let SUF], FORMATOS);
end
end

% =========================================================================
% Convergencia
% =========================================================================
if hacer(4)
printf('=== 4. Convergencia ===\n');  fflush(stdout);

Eu = zeros(numel(Ns), 2);  Ea = zeros(numel(Ns), 2);
for s = 1:2
  p = P(s);
  for m = 1:numel(Ns)
    xu = linspace(0, 1, Ns(m)+1);
    Uu = resuelveRobin1D(xu, p.k, p.f, p.ga, p.gb);
    Eu(m,s) = erroresMalla(xu, Uu, p.ue);
    [xa, Ua] = adaptaDeBoor(xu, p.k, p.f, p.ga, p.gb, oA);
    Ea(m,s) = erroresMalla(xa, Ua, p.ue);
    printf('  %s  N=%4d   uniforme %.4e   adaptada %.4e   ganancia %.1f\n', ...
           p.let, Ns(m), Eu(m,s), Ea(m,s), Eu(m,s)/Ea(m,s));
    fflush(stdout);
  end
end

for s = 1:2
  f = nuevaFigura(ANCHO, ALTO, LETRA);
  loglog(Ns, Eu(:,s), 'ks-', 'linewidth', 1.2, 'markersize', 5); hold on;
  loglog(Ns, Ea(:,s), 'bo-', 'linewidth', 1.4, 'markersize', 5);
  loglog(Ns, Ea(1,s)*(Ns(1)./Ns).^2, 'k:', 'linewidth', 1.2);
  xlabel('N'); ylabel(R.err); grid on;
  leyenda(LETRA-1, 'southwest', R.uniforme, R.adaptada, R.pend2);
  guarda(f, dirFig, ['fig_convergencia_' P(s).let SUF], FORMATOS);
end
end

printf('\nFiguras escritas en %s\n', dirFig);
end

% ==========================================================================
function f = nuevaFigura(w, h, fs)
% Figura con el motor, el tamaño de papel y el cuerpo de letra fijados. El
% papel ha de coincidir con la posición para que el PDF salga recortado al
% dibujo y sin márgenes.
  if nargin < 3 || isempty(fs), fs = 8; end
  f = figure('visible', 'off', 'position', [0 0 w h]);
  if any(strcmp(available_graphics_toolkits(), 'gnuplot'))
    set(f, '__graphics_toolkit__', 'gnuplot');
  end
  set(f, 'paperunits', 'points', 'papersize', [w h], ...
         'paperposition', [0 0 w h]);
  set(f, 'defaultaxesfontsize', fs, 'defaulttextfontsize', fs, ...
         'defaultlinelinewidth', 1);
end

function h = leyenda(fs, loc, varargin)
% Leyenda con el cuerpo de letra fijado a mano. Hace falta porque el motor de
% gnuplot no hereda defaulttextfontsize en la leyenda, y a este tamaño de panel
% la caja se come un tercio del dibujo. 'best' no está implementado en gnuplot,
% así que la esquina se elige aquí.
  h = legend(varargin{:}, 'location', loc);
  set(h, 'fontsize', fs);
end

function guarda(f, dirFig, nombre, formatos)
% Escribe la figura en todos los formatos pedidos y la cierra.
  for m = 1:numel(formatos)
    d = formatos{m};
    switch d
      case '-dpng',  print(f, fullfile(dirFig, [nombre '.png']), d, '-r140');
      case '-dpdf',  print(f, fullfile(dirFig, [nombre '.pdf']), d);
      case '-dsvg',  print(f, fullfile(dirFig, [nombre '.svg']), d);
      case '-depsc2', print(f, fullfile(dirFig, [nombre '.eps']), d);
      otherwise, warning('guarda: formato %s desconocido', d);
    end
  end
  close(f);
end

function c = centros(x)
  c = (x(1:end-1) + x(2:end))/2;
end

function s = salto(x)
  h = diff(x);
  s = max(abs(h(2:end)./h(1:end-1) - 1));
end

% Los motores qt y fltk rasterizan con OpenGL y exigen ventana visible, así
% que fallan con octave-cli sin servidor gráfico, y tampoco escriben vectorial.
% gnuplot escribe PDF, PNG y SVG directamente. Hay que llamar a esto antes de
% crear las figuras, porque graphics_toolkit sólo afecta a las figuras nuevas;
% aun así nuevaFigura vuelve a fijarlo en cada una, porque el ajuste global no
% siempre alcanza a las que se crean después.
function tk = backendGrafico()
  tmp = [tempname() '.png'];
  tk  = graphics_toolkit();
  ok  = probar(tmp);
  if ~ok && any(strcmp(available_graphics_toolkits(), 'gnuplot'))
    warning('off', 'Octave:gnuplot-graphics');
    graphics_toolkit('gnuplot');
    tk = 'gnuplot';
    ok = probar(tmp);
  end
  if exist(tmp, 'file'), delete(tmp); end
  if ~ok
    warning('backendGrafico: ningún motor gráfico puede escribir PNG');
  end
end

function ok = probar(archivo)
  ok = true;
  f = figure('visible', 'off');
  try
    plot(1:3, 1:3);
    print(f, archivo, '-dpng');
  catch
    ok = false;
  end
  close(f);
end
