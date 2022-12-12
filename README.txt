El programa 'DBCleaner' sirve para filtrar las bases de datos según la información que aparece en la descripcion del fasta de las bases de datos.
También elimina las secuencias que sean exactamente iguales (tanto en secuencia como en longitud). En el caso de que dos o más secuencias sean iguales
la descripción del fasta resultante será la unión de los códigos Accession de cada uno separado por '__'.

Para ejecutar el programa se requiere de un archivo params.yalm donde se indica la ruta al archivo de la base de datos, y los filtros que se le va a aplicar.
Hay dos tipos de filtros:
	-OR: si se indica más de una palabra buscará que aparecezca una o las otras.
	-AND: deberá aparecen todas las palabras indicadas en la descripción
Si no se quiere aplicar un filtro se tiene que poner None.
También se indica en el params el nombre del archivo resultante.
Un ejemplo de archivo params sería:

############################

fasta_1:
  - file: path
  - filter_OR:
    - None
  - filter_AND:
    - protein_coding
    - gene_symbol


fasta_2:
  - file: path
  - filter_OR:
    - PE=1
    - PE=2
  - filter_AND:
    - None

fasta_3:
  - file: path
  - filter_OR:
    - NP_
  - filter_AND:
    - None

outfile: path

#############################

El programa se ejecuta desde la cmd poniendo:

>python path_to_program -p path_to_params