% prueba_truncamiento.m
% Comprueba los tres apartados de la proposicion sobre el truncamiento del
% suavizado eliptico, y la franja que queda entre la positividad de los pesos
% y la convergencia de la serie.
%
% La version discreta de (I - beta_s^{-2} d^2/dxi^2) rho~ = rho sobre malla
% computacional uniforme de paso dxi es
%
%     (I - t T) rho~ = rho,      t = 1/(beta_s^2 dxi^2),
%
% con (T v)_i = v_{i-1} - 2 v_i + v_{i+1}. El promedio de tres puntos que se
% obtiene al truncar la serie de Neumann en el primer termino es I + t T.
%
% Lo que se verifica:
%
%   (i)   los pesos (t, 1-2t, t) son no negativos si y solo si t <= 1/2,
%         es decir beta_s^2 dxi^2 >= 2;
%   (ii)  el radio espectral de t T tiende a 4t, luego la serie converge
%         uniformemente en el numero de celdas si y solo si 4t < 1, es decir
%         beta_s^2 dxi^2 > 4;
%   (iii) sobre el modo de signos alternos v_i = (-1)^i, que en el interior
%         cumple T v = -4 v, el promedio truncado multiplica por 1 - 4t y el
%         operador sin truncar por 1/(1 + 4t).
%
% La franja 2 <= beta_s^2 dxi^2 < 4 es la interesante: los pesos son
% admisibles y la serie ya diverge, de modo que el promedio le invierte el
% signo al modo de signos alternos en vez de amortiguarlo. En el extremo
% beta_s^2 dxi^2 = 2 se lo devuelve con la amplitud intacta.
%
% Las sumas parciales que imprime la ultima tabla son las que se citan en el
% articulo: -1, 3, -5, 11, -21, 43, -85.

function prueba_truncamiento()

n = 400;
e = ones(n, 1);
T = spdiags([e, -2*e, e], -1:1, n, n);
v = (-1).^(1:n).';          % modo de signos alternos
c = round(n/2);             % indice interior, lejos de los dos bordes

printf('\n=== espectro de T ===\n');
lam = eig(full(T));
printf('  autovalores en [%.6f, %.6f]      teoria: (-4, 0)\n', ...
       min(lam), max(lam));
printf('  (T v)/v en el interior = %.6f      teoria: -4\n', (T*v)(c)/v(c));

printf('\n=== los tres apartados, valor por valor ===\n');
printf('  %-8s %-8s %-16s %-5s %-8s %-6s %-9s %-9s\n', ...
       'bs2dxi2', 't', 'pesos', 'pos', 'radio', 'serie', '3 puntos', 'exacto');
for bd2 = [1.5 2 2.5 3 4 5 6 10]
  t  = 1/bd2;
  ap = (speye(n) + t*T)*v;          % promedio de tres puntos
  ex = (speye(n) - t*T)\v;          % operador sin truncar
  if 4*t < 1, s = 'conv'; else, s = 'div'; end
  printf('  %-8.1f %-8.4f (%+.3f,%+.3f) %-5d %-8.4f %-6s %+-9.4f %+-9.4f\n', ...
         bd2, t, t, 1-2*t, 1-2*t >= 0, 4*t, s, ap(c)/v(c), ex(c)/v(c));
end
printf('\n  La columna «pos» marca la positividad de los pesos y «serie» la\n');
printf('  convergencia. Difieren en la franja 2 <= bs2dxi2 < 4, donde el\n');
printf('  promedio de tres puntos da un factor negativo.\n');

printf('\n=== sumas parciales de la serie en bs2dxi2 = 2 ===\n');
t = 0.5;  s = v;  acc = v;
printf('  coeficiente del modo de signos alternos:');
for m = 1:7
  s   = t*(T*s);
  acc = acc + s;
  printf(' %g', acc(c)/v(c));
end
printf('\n  El operador sin truncar da 1/(1+4t) = %.4f.\n\n', 1/(1+4*t));

end
