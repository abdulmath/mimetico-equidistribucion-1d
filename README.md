# Adaptación de malla por equidistribución para diferencias finitas miméticas

Código en GNU Octave que acompaña al artículo *La escala del suavizado en la
adaptación de malla por equidistribución para diferencias finitas miméticas*,
de Abdul Abner Lugo Jiménez (IEES Loyola, San Cristóbal, República Dominicana).

El repositorio contiene el esquema mimético de Castillo y Grone de orden dos en
una dimensión, el lazo de adaptación por equidistribución con el algoritmo de de
Boor, los tres suavizados de la función monitor que se comparan en el artículo,
y los programas que generan todas las tablas y figuras del texto. Ningún número
del artículo está copiado a mano: todos salen de estos programas.

## Requisitos

GNU Octave 8 o posterior. No hace falta ningún paquete adicional. Las figuras se
escriben con el motor gnuplot, que Octave ya trae, para poder generarlas sin
servidor gráfico.

## Cómo reproducir los resultados

Cada programa se ejecuta tal cual, por ejemplo

```sh
cd pruebas
octave-cli --eval "prueba_dos_problemas"
```

Para rehacer de una vez todos los archivos de `salidas/`:

```sh
sh regenera_salidas.sh
```

La única sensibilidad apreciable de las cifras es la tolerancia con que se
detiene la iteración exterior. En el problema A con 3200 celdas el error en
norma del máximo vale 7.7364e-07 con `tolX = 1e-10` y 7.7366e-07 con `1e-9`,
porque a esa altura la malla todavía se mueve en la novena cifra. Todos los
programas usan 1e-10 salvo `prueba_suavizado_global.m`, que usa 1e-9 por el
motivo que se explica más abajo.

El número de hilos de la biblioteca de álgebra lineal no cambia ninguna de las
cifras publicadas; se comprobó ejecutando `prueba_dos_problemas` con los hilos
fijados a uno y con los hilos libres, y las dos corridas coinciden dígito a
dígito. `regenera_salidas.sh` los fija de todos modos, para que cualquier
diferencia que aparezca sea fácil de localizar.

Cada programa de `pruebas/` es independiente y añade `codigo/` a la ruta por su
cuenta, así que puede ejecutarse desde cualquier sitio mientras `codigo/` y
`pruebas/` sigan siendo carpetas hermanas.

La carpeta `salidas/` guarda la salida de texto de cada programa tal como se
generó para el artículo. Sirve para comparar sin tener que volver a correr las
mallas grandes, que tardan.

## Qué genera cada programa

| Programa | Qué produce |
|---|---|
| `prueba_identidad.m` | Tabla 1: recuperación de P en la identidad mimética sobre mallas no uniformes |
| `prueba_dos_problemas.m` | Tablas 2, 3 y 4: malla uniforme y malla adaptada en los problemas A y B |
| `prueba_suavizado_global.m` | Tablas 5 y 7: comparación de los tres suavizados y barrido de beta_s |
| `prueba_limite_filtro.m` | Tabla 6: límite de la densidad suavizada por dos caminos independientes |
| `prueba_truncamiento.m` | Proposición sobre el truncamiento del suavizado elíptico |
| `prueba_lazo_interno.m` | Que la recurrencia h_i = C/rho_i no converge a la malla equidistribuida |
| `prueba_deboor.m` | Que un paso de de Boor equidistribuye hasta el redondeo; dependencia del factor de relajación |
| `prueba_equidistribucion_nueva.m` | Convergencia del lazo exterior, movimiento máximo iteración a iteración |
| `prueba_relajacion.m` | Efecto del factor de subrelajación en mallas gruesas, en los dos problemas |
| `prueba_ganancia.m` | Una etapa más allá de las tablas, hasta 6400 celdas: la ganancia se estabiliza |
| `prueba_tolerancia.m` | Por qué la comparación del suavizado elíptico se detiene con 1e-9 y no con 1e-10 |
| `figuras_articulo.m` | Los diez paneles de las figuras, en PDF y PNG, dentro de `figuras/` |

Los programas `prueba_orden.m` a `prueba_orden7.m` son la secuencia de
diagnóstico que llevó al criterio sobre la escala del suavizado, en el orden en
que se hicieron. `prueba_orden.m` y `prueba_orden3.m` descartan las dos
explicaciones alternativas que discute el artículo, la transformación de la
variable y el tratamiento de la frontera. `prueba_orden2.m` es el experimento
de control, con la densidad construida al margen de la malla a partir de la
derivada exacta. `prueba_orden4.m` compara el filtro de anchura fija en x con
el de Dorfi y Drury y barre el parámetro de escala del monitor.
`prueba_orden5.m` barre sigma_0 y `prueba_orden6.m` y `prueba_orden7.m`
estudian la anchura adaptativa y la relajación del lazo en mallas gruesas.

`figuras_articulo.m` admite un argumento con las etapas que se quieren rehacer,
por ejemplo `figuras_articulo([3 4])` para el salto y la convergencia.

## Qué hay en `codigo/`

| Archivo | Función |
|---|---|
| `operadoresMimeticos1D.m` | Gradiente, divergencia y pesos de Castillo y Grone sobre malla mimética no uniforme |
| `resuelveRobin1D.m` | Resuelve el problema de contorno con condiciones tipo Robin |
| `monitorArco.m` | Función de densidad de malla de longitud de arco escalada |
| `suavizaAnchuraFija.m` | Promedio gaussiano de escala fija en x |
| `suavizaDorfiDrury.m` | Promedio de escala fija en número de celdas |
| `suavizaGlobal.m` | Suavizado elíptico de Huang y Russell |
| `deBoor.m` | Inversión exacta de la primitiva lineal a trozos, O(N) |
| `equidistribuye.m` | Un paso de equidistribución |
| `adaptaDeBoor.m` | Lazo exterior de adaptación con subrelajación |
| `adaptaEquidistribucion.m` | Variante del lazo con la recurrencia clásica, para comparar |
| `erroresMalla.m` | Errores en norma del máximo y en norma discreta L2 |

## Advertencia sobre la tolerancia

Las corridas del suavizado elíptico se detienen con una tolerancia de 1e-9, un
orden por encima de la que usa el resto del artículo. Con 1e-10 el problema A
con 3200 celdas y beta_s = 20 entra en ciclo. El error final no cambia en
ninguna de las cifras publicadas, pero la malla queda determinada un orden de
magnitud peor que con el filtro de escala fija, que sí alcanza 1e-10.

## Cita

Si usa este código, cite el artículo. El archivo `CITATION.cff` tiene los datos
en forma legible por máquina.

## Licencia

MIT, en `LICENSE`.
