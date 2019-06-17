clear
name1=datos
name2=figuras

mkdir "$name1" "$name2"
gcc -Wall -O2 Hartmann2.cpp -o H
 
for k in {-1e-3,-1e-2,-1e-1,-1}
do
 echo "Se está calculando el campo magnético B = $k"
 ./H "$k"
 echo "$k" | python3 CambioNombre.py
 mv *"$k".dat "${name1}"
done

