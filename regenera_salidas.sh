#!/bin/sh
# regenera_salidas.sh
# Regenera todos los archivos de salidas/ con la biblioteca de álgebra lineal
# en un solo hilo, que es la condición en la que se obtuvieron las cifras
# publicadas.
#
#   sh regenera_salidas.sh            todos los programas
#   sh regenera_salidas.sh prueba_X   sólo ese
#
# El producto denso que evalúa el filtro de escala fija se suma en un orden
# que depende del número de hilos, y ese orden cambia la quinta cifra
# significativa de los errores. OMP_NUM_THREADS no basta por sí solo: OpenBLAS
# atiende antes a OPENBLAS_NUM_THREADS, y las compilaciones contra MKL, a
# MKL_NUM_THREADS. Se fijan las tres.

set -e

OMP_NUM_THREADS=1
OPENBLAS_NUM_THREADS=1
MKL_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1
export OMP_NUM_THREADS OPENBLAS_NUM_THREADS MKL_NUM_THREADS VECLIB_MAXIMUM_THREADS

raiz=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$raiz/salidas"

if [ $# -gt 0 ]; then
  lista="$*"
else
  lista="prueba_identidad
         prueba_dos_problemas
         prueba_suavizado_global
         prueba_limite_filtro
         prueba_truncamiento
         prueba_lazo_interno
         prueba_deboor
         prueba_equidistribucion_nueva
         prueba_relajacion
         prueba_ganancia
         prueba_tolerancia
         prueba_orden
         prueba_orden2
         prueba_orden3
         prueba_orden4
         prueba_orden5
         prueba_orden6
         prueba_orden7"
fi

for p in $lista; do
  printf '%s ... ' "$p"
  if octave-cli --quiet --path "$raiz/pruebas" --eval "$p" \
       > "$raiz/salidas/$p.txt" 2>"$raiz/salidas/$p.err"; then
    if [ -s "$raiz/salidas/$p.err" ]; then
      echo "hecho, con avisos en salidas/$p.err"
    else
      rm -f "$raiz/salidas/$p.err"
      echo "hecho"
    fi
  else
    echo "FALLÓ, ver salidas/$p.err"
  fi
done

echo
echo "Comprobación rápida de la cifra sensible:"
grep -h "3200" "$raiz/salidas/prueba_dos_problemas.txt" 2>/dev/null | head -3
echo "(el problema A con 3200 celdas y malla adaptada debe dar 7.7364e-07)"
