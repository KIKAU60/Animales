import streamlit as st
from Bio import SeqIO
import matplotlib.pyplot as plt

# Función para leer un archivo FASTA y obtener la secuencia
def read_fasta(file):
    try:
        # Leer el archivo FASTA
        record = SeqIO.read(file, "fasta")
        sequence = record.seq
        description = record.description
        return sequence, description
    except Exception as e:
        return None, f"Error al leer el archivo FASTA: {e}"

# Función para calcular el contenido GC de la secuencia
def calculate_gc_content(sequence):
    if sequence:
        sequence = str(sequence).upper()  # Convertir la secuencia a string y en mayúsculas
        gc_count = sum(1 for base in sequence if base in 'GC')
        return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0
    return 0  # Si la secuencia está vacía o no definida, devolver 0

# Función para visualizar la secuencia y otras informaciones
def visualize_gene_info(sequence, description):
    if sequence and len(sequence) > 0:
        sequence = str(sequence).upper()  # Convertir la secuencia a mayúsculas
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)

        # Mostrar la descripción del gen
        st.subheader("Descripción del Gen:")
        st.write(description)

        # Mostrar información adicional
        st.subheader("Información sobre la Secuencia:")
        st.write(f"Longitud de la secuencia: {length} bases")
        st.write(f"Contenido GC: {gc_content:.2f}%")
        
        # Mostrar gráfico de contenido GC
        st.subheader("Gráfico de distribución de bases")
        fig, ax = plt.subplots()
        ax.bar(["A", "T", "G", "C"], [sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C')])
        ax.set_ylabel("Número de bases")
        ax.set_title("Distribución de bases en la secuencia de ADN")
        st.pyplot(fig)

    else:
        st.warning("La secuencia de ADN no se pudo obtener o está vacía. Por favor, verifica el archivo FASTA.")

# Título de la aplicación
st.title("Análisis de Secuencia de ADN desde un archivo FASTA")

# Cargar archivo FASTA
uploaded_file = st.file_uploader("Sube tu archivo FASTA", type=["fasta"])

# Botón para procesar el archivo
if uploaded_file is not None:
    sequence, description = read_fasta(uploaded_file)
    if sequence:
        st.success("Secuencia de ADN cargada exitosamente!")
        visualize_gene_info(sequence, description)
    else:
        st.error(description)
else:
    st.warning("Por favor, sube un archivo FASTA para continuar.")
