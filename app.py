import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
import seaborn as sns
import random
import re
from Bio import Entrez
from Bio import SeqIO
from Bio import Phylo
import nglview as nv

# Configuración de Entrez: Debes proporcionar tu correo para cumplir con los requisitos de GenBank
Entrez.email = "a223201128@unison.mx"  # Sustituye con tu correo electrónico

# Función para obtener la secuencia de ADN desde GenBank
def obtener_secuencia_genbank(accession_number):
    try:
        # Fetch la secuencia de GenBank usando el número de acceso
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")  # Lee el archivo de GenBank
        secuencia_adn = str(record.seq)  # Convierte la secuencia en cadena
        return secuencia_adn, record
    except Exception as e:
        st.error(f"Error al obtener la secuencia de GenBank: {e}")
        return None, None

# Función para calcular las proporciones de nucleótidos (A, T, C, G)
def calcular_proporcion_nucleotidos(secuencia):
    counter = Counter(secuencia)
    total = len(secuencia)
    proporciones = [counter['A'] / total, counter['T'] / total, counter['C'] / total, counter['G'] / total]
    return proporciones

# Función para graficar la frecuencia de codones (tripletas) usando Plotly
def graficar_codones_interactivo(secuencia):
    codones = [secuencia[i:i+3] for i in range(0, len(secuencia), 3)]
    counter = Counter(codones)
    codones_unicos = list(counter.keys())
    frecuencias = list(counter.values())

    # Gráfico interactivo usando Plotly
    fig = go.Figure([go.Bar(x=codones_unicos, y=frecuencias, marker_color='skyblue')])
    fig.update_layout(
        title="Frecuencia de Codones en la Secuencia de ADN",
        xaxis_title="Codón",
        yaxis_title="Frecuencia",
        xaxis_tickangle=-45
    )
    st.plotly_chart(fig)

# Función para generar un árbol filogenético de secuencias de ADN
def generar_arbol_filogenetico(secuencias):
    # Cargar las secuencias de ADN y construir un árbol filogenético
    from Bio.Align import MultipleSeqAlignment
    from Bio.Align.Applications import ClustalwCommandline

    # Crear un alineamiento ficticio de ejemplo con las secuencias proporcionadas
    alignment = MultipleSeqAlignment([SeqIO.read(seq, "genbank") for seq in secuencias])

    # Usar ClustalW para generar un árbol filogenético
    clustalw_cline = ClustalwCommandline("clustalw2", infile="alignment.aln")
    clustalw_cline()

    # Leer el archivo de salida y graficar el árbol filogenético
    tree = Phylo.read("tree.dnd", "newick")
    Phylo.draw(tree)
    st.pyplot()

# Función para mostrar la estructura 3D de la hélice de ADN
def mostrar_helice_3d(secuencia):
    # Visualización 3D de la hélice de ADN usando nglview
    estructura = f'PDB_file_{random.randint(1, 10000)}.pdb'  # Generamos un archivo temporal para la hélice
    # Aquí deberías usar una estructura de ADN en formato PDB, por ahora usaremos un nombre genérico

    # Muestra la estructura 3D
    view = nv.show_file(estructura)  # Asume que tienes un archivo PDB para la hélice
    view

# Función principal para la aplicación Streamlit
def main():
    st.title("Análisis de Secuencias de ADN desde GenBank")
    
    # Entrada para el número de acceso de GenBank (Accession Number)
    accession_number = st.text_input("Introduce el número de acceso de GenBank (Accession Number):")
    
    if accession_number:
        # Obtener la secuencia de ADN desde GenBank
        secuencia_adn, record = obtener_secuencia_genbank(accession_number)
        
        if secuencia_adn:
            st.write(f"Secuencia de ADN obtenida de GenBank (Accession Number: {accession_number}):")
            st.text(secuencia_adn)
            
            # Opciones de ilustración
            ilustracion = st.selectbox(
                "Selecciona la ilustración para la secuencia de ADN:",
                ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Árbol Filogenético', 'Estructura 3D de la Hélice']
            )
            
            # Actualizar visualización según la opción seleccionada
            if ilustracion == 'Proporciones de Nucleótidos':
                proporciones = calcular_proporcion_nucleotidos(secuencia_adn)
                fig, ax = plt.subplots(figsize=(8, 8))
                ax.pie(proporciones, labels=['A', 'T', 'C', 'G'], autopct='%1.1f%%', startangle=140)
                ax.set_title("Proporción de Nucleótidos en la Secuencia de ADN")
                st.pyplot(fig)

            elif ilustracion == 'Frecuencia de Codones':
                graficar_codones_interactivo(secuencia_adn)

            elif ilustracion == 'Árbol Filogenético':
                # Asegúrate de tener las secuencias necesarias para el árbol filogenético
                secuencias = [secuencia_adn]  # Aquí debes añadir más secuencias si las tienes
                generar_arbol_filogenetico(secuencias)

            elif ilustracion == 'Estructura 3D de la Hélice':
                mostrar_helice_3d(secuencia_adn)
        else:
            st.write("No se pudo obtener la secuencia. Verifica el número de acceso.")
    
# Ejecutar la aplicación
if __name__ == "__main__":
    main()
