% prueba_relajacion.m
% Efecto del factor de subrelajación sobre la convergencia del lazo exterior
% en mallas gruesas, con el filtro de anchura fija en x.
%
% Con omega = 1 el lazo converge mientras la malla resuelva la capa. Cuando no
% la resuelve entra en un ciclo y agota el tope de iteraciones. La relajación
% lo saca del ciclo sin cambiar el punto fijo: para cada tamaño de malla los
% valores de omega que convergen dan el mismo error en las cinco cifras que se
% reportan, y sólo cambia el número de iteraciones.
%
% Genera las cifras de la sección "La iteración exterior y su subrelajación".

function prueba_relajacion()
raiz = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(raiz, 'codigo'));

% --- Problema A : pico interior, coeficiente variable ---------------------
al = 250;  x0 = 0.8;
A.et = 'A';
A.k  = @(z) 1/al + al*(z-x0).^2;
A.f  = @(z) -2*(1 + al*(z-x0).*(atan(al*x0) + atan(al*(z-x0))));
A.ue = @(z) (1-z).*(atan(al*(z-x0)) + atan(al*x0));
A.ga = al/(1+al^2*x0^2);
A.gb = -atan(al*(1-x0)) - atan(al*x0);

% --- Problema B : capa en el borde izquierdo, coeficiente constante -------
c = atan(100);
B.et = 'B';
B.k  = @(z) ones(size(z));
B.f  = @(z) -2e6*z./(c*(1 + 1e4*z.^2).^2);
B.ue = @(z) atan(100*z)/c;
B.ga = 100/c;
B.gb = 1 + 100/(c*(1+1e4));

printf('\n=== Relajacion en mallas gruesas, filtro de anchura fija ===\n');
printf('  %-4s %-8s %-7s %-8s %-8s %-13s\n', ...
       'Prob', 'Celdas', 'omega', 'conv', 'iter', 'Linf');

for P = [A B]
  for N = [50 100 200]
    x0m = linspace(0, 1, N+1);
    for om = [1.00 0.90 0.70 0.50]
      o = struct('omega', om, 'maxIt', 400, 'tolX', 1e-10);
      [xa, Ua, h] = adaptaDeBoor(x0m, P.k, P.f, P.ga, P.gb, o);
      eI = erroresMalla(xa, Ua, P.ue);
      if h.convergio, s = 'si'; else, s = 'NO'; end
      printf('  %-4s %-8d %-7.2f %-8s %-8d %-13.4e\n', ...
             P.et, N, om, s, h.it, eI);
    end
  end
end
printf('\n');
end
