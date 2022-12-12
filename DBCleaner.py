######### Module metadata variables ########
__author__ = "David del Rio Aledo"
__credits__ = ["David del Rio", "Jose Manuel Rodriguez", "Inmaculada Jorge"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "David del Rio"
__status__ = "Development"

###########################################



                #########################
                # Import global modules #
                #########################


#################  Tkinter libraries  #################
from tkinter import *
from tkinter import ttk
######################################################


#################  Others libraries  #################
from Bio import SeqIO
import pandas as pd
import copy
import yaml
from yaml import Loader, Dumper
import argparse
from progress.bar import Bar, ChargingBar
######################################################


####################
# Global variables #
####################

sequences = []
descriptions=[]
ids=[]

fastas=[]

sequences_filtered=[]
sequences_exclude=[]

'''
list_entries=[]
list_labels_files=[]
list_labels_columns=[]


#### Table entries size ###
total_rows=3
total_columns=2

### Output information ###
output_iformation="Datos del procesado de la base de datos \n \n"
'''


####################
####################




###################
# Root activation #
###################
'''
root= Tk()
root.geometry("700x350")
'''




######################
# Output files names #
######################

#output_file="C:/Users/ddelrioa/Desktop/Database_PE1_tr_y_sw/prueba_2.fasta"
#output_file_100="C:/Users/ddelrioa/Desktop/Database_PE1_tr_y_sw/prueba_2_100.fasta"



#####################
# Filters functions #
#####################

def clean_fasta(df,sequences,ids,descriptions):
    #For each file check if it has the text put in filters of params file.
    #There are two types of filters: OR or AND, if you don't want someone put None in filters params.
    #global output_iformation

    index_df=-1

    #with open(output_file,"w") as out_file:

    for file in df[0]:
        print(file)
        count_total=0
        count_filtred=0

        index_df+=1
        arguments_OR=df[1][index_df].split(sep=',')
        arguments_AND=df[2][index_df].split(sep=',')

        for fasta in SeqIO.parse(file,"fasta"):

            name, sequence, description=fasta.id, str(fasta.seq), fasta.description
            count_total+=1
            filtered="No"

            for flt in arguments_OR:

                if str(flt) != None and str(flt) in description:

                    #out_file.write(">" + name + "\n" + sequence + "\n")
                    sequences.append(sequence)
                    ids.append(name)
                    descriptions.append(description)
                    fastas.append(fasta)
                    count_filtred+=1
                    filtered="Yes"

            if arguments_AND !=None and all(x in description for x in arguments_AND)==True:

                    #out_file.write(">" + name + "\n" + sequence + "\n")
                    sequences.append(sequence)
                    ids.append(name)
                    descriptions.append(description)
                    fastas.append(fasta)
                    count_filtred+=1
                    filtered="Yes"


            if filtered=="No":
                sequences_exclude.append(name)

        print("El archivo " +str(index_df+1)+ " tiene " + str(count_total) +" secuencias totales")
        print("El archivo " +str(index_df+1)+ " tiene " + str(count_filtred) +" secuencias filtradas\n")

    print("En total se han eliminado " + str(len(sequences_exclude)) +" secuencias\n\n")

    '''
        output_iformation=output_iformation+"El archivo " +str(index_df)+ " tiene " + str(count_total) +" secuencias totales\n"
        output_iformation=output_iformation + "El archivo " +str(index_df)+ " tiene " + str(count_filtred) +" secuencias filtradas\n"
        change_information(text=output_iformation)
        
    output_iformation=output_iformation+ "En total se han eliminado " + str(len(sequences_exclude)) +" secuencias\n\n"
    change_information(text=output_iformation)
    '''
    

def filterer(fastas,descriptions,outfile):
    #Check if there are equal sequences after applying clean_fasta function (fastas object is filled in that function)
    #Those that are the same unite the accesions and attach all unique sequences to the file.
    #global output_iformation

    count_total=0
    count_equal=0
    fastas_comparation=copy.copy(fastas)
    bar=ChargingBar('Comparando secuencias y eliminando repetidas', max=int(len(fastas)))

    with open(outfile,"w") as out_file_100:

        for fasta in fastas:

            name, sequence, description=fasta.id, str(fasta.seq), fasta.description
            fastas_comparation.pop(0)
            count_total+=1

            name_id=copy.copy(name)

            equal=0

            for i in range(len(fastas_comparation)):

                if sequence == str(fastas_comparation[i].seq):

                    equal=1

                    sequences_filtered.append(sequence)
                    name_id=name_id+"__"+str(fastas_comparation[i].id)
                    count_equal+=1
                    #name_accession_1=re.search(r'[\|]([\w]+[\|])',name)

            if equal==0:
                out_file_100.write(">" + name + "\n" + sequence + "\n")

            else:
                out_file_100.write(">" + name_id + "\n" + sequence + "\n")

            bar.next()

        bar.finish()
        print("De un total de " +str(count_total)+ " secuencias, son iguales " + str(count_equal)+"\n")
        #output_iformation=output_iformation+"De un total de " +str(count_total)+ " secuencias, son iguales " + str(count_equal)+"\n"
        #change_information(text=output_iformation)


#def change_information(text):
#    otput_information_lbl.configure(text=text)

def params_reader(params):
    #Read the params input file and transform to dataframe.
    #Input params files, return dataframe and path of outfile

    with open(params,'r') as params:
        data=yaml.load(params,Loader=Loader)

    frame_data=[]
    index_frame_data=-1

    for clave in data:
        index_frame_data+=1

        if clave!='outfile':
            frame_data.append([])
            frame_data[index_frame_data].append(data[clave][0]['file'])
            frame_data[index_frame_data].append(','.join(data[clave][1]['filter_OR']))
            frame_data[index_frame_data].append(','.join(data[clave][2]['filter_AND']))

        else:
            outfile_path=data[clave]

    df=pd.DataFrame(frame_data)
    print(df)
    print('\n')

    return df, outfile_path


#################
# Main function #
#################

def main(args):
    #Then execute clean_fasta and filterer functions.

    df, outfile_path=params_reader(params=args.paramsfile)
    clean_fasta(df,sequences,ids,descriptions)

    filterer(fastas,descriptions,outfile_path)



##########
# Parser #
##########


parser=argparse.ArgumentParser(prog='DBCleaner',
    description='\n\n\nFilter fatas in each input file with the filters indicated searching the text in fasta description.\n'+ 
    'Then compare all fasta filtreted to join the same sequences as one, the final accession is the union of those who are equal\n\n\n'
    )
parser.add_argument('-p','--paramsfile',required=True, help='Input params file with all path to database files and filters')
args=parser.parse_args()


main(args)



'''
def main(list_entries):
    #Read the text inserted in entries and create a dataframe.
    #Then execute clean_fasta and filterer functions.

    frame_data=[]

    for i in range(total_rows):
        frame_data.append([])

        for j in range(total_columns):
            frame_data[i].append(list_entries[i][j].get())

    df=pd.DataFrame(frame_data)

    print(df)
    clean_fasta(df,sequences,ids,descriptions)

    filterer(fastas,df,descriptions)


################
# GUI creation #
################


###    Entries and labels   ###

for i in range(total_rows):

    File_lbl=Label(root,text="File " + str(i+1))
    File_lbl.configure(bg="white", width=10, relief="groove", bd=2, fg="#000000")
    File_lbl.grid(row=i+1, column=0)
    list_labels_files.append(File_lbl)

for j in range(total_columns):

    if j == 0:
        column_lbl=Label(root,text="File path")
        column_lbl.configure(bg="white", width=17, relief="groove", bd=2, fg="#000000")
        column_lbl.grid(row=0, column=j+1)
        list_labels_columns.append(column_lbl)

    else:
        column_lbl=Label(root,text="Filter " + str(j))
        column_lbl.configure(bg="white", width=17, relief="groove", bd=2, fg="#000000")
        column_lbl.grid(row=0, column=j+1)
        list_labels_columns.append(column_lbl)

for i in range(total_rows):
    list_entries.append([])

    for j in range(total_columns):
        entry=Entry(root, width=20, fg="black",bg="white")
        entry.grid(row=i+1,column=j+1)
        list_entries[i].append(entry)


otput_information_lbl=Label(root,text=output_iformation)
otput_information_lbl.configure(bg="white", width=50, relief="groove", bd=2, fg="#000000")
otput_information_lbl.place(x=500, y=10)

        
###    Buttons   ###

start_btn=Button(root, text="Start",command=lambda:[main(list_entries=list_entries)],
    bg="grey",state="normal")
start_btn.place(x=400,y=30)



root.mainloop()
'''