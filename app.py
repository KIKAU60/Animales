import streamlit as st
from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt

# Configuración de Entrez
Entrez.email = "tu_email@example.com"  # Asegúrate de poner tu correo electrónico

# Función para obtener la secuencia de ADN desde GenBank usando un Accession ID
def fetch_dna_from_genbank(accession_id):
    try:
        # Buscar en GenBank utilizando el Accession ID
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Verificar si la secuencia es válida
        if hasattr(record, 'seq') and record.seq:
            return record.seq, None
        else:
            return None, "Secuencia no encontrada en el archivo de GenBank."

    except Exception as e:
        return None, f"Error al obtener la secuencia de GenBank: {e}"

# Función para calcular el contenido GC de la secuencia
def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in 'GCgc')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

# Función para analizar y mostrar la secuencia de ADN
def visualize_dna(sequence):
    if sequence:
        sequence = sequence.strip().upper()  # Limpiar y convertir a mayúsculas
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)

        # Mostrar la secuencia de ADN
        st.subheader("Secuencia de ADN:")
        st.text(sequence)

        # Mostrar información adicional
        st.subheader("Información adicional:")
        st.write(f"Longitud de la secuencia: {length} bases")
        st.write(f"Contenido GC: {gc_content:.2f}%")
        st.write("Posibles características genéticas:")
        st.write("- Aumento en la producción de proteínas (si es un gen codificador).")
        st.write("- Variante genética de adaptación a diferentes ambientes.")

        # Mostrar gráfico de contenido GC
        st.subheader("Gráfico del contenido GC")
        fig, ax = plt.subplots()
        ax.bar(["A", "T", "G", "C"], [sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C')])
        ax.set_ylabel("Número de bases")
        ax.set_title("Distribución de bases en la secuencia de ADN")
        st.pyplot(fig)

    else:
        st.warning("Por favor, introduce una secuencia de ADN válida.")

# Título de la aplicación
st.title("Análisis de Secuencias de ADN Animal")

# Entrada para el ID de GenBank
accession_id = st.text_input("Introduce el ID de acceso de GenBank (Accession ID):", 
                             help="Puedes buscar el ID de acceso de un gen o genoma en GenBank (por ejemplo, NC_000913)")

# Botón para obtener la secuencia desde GenBank
if st.button("Obtener secuencia desde GenBank"):
    if accession_id:
        sequence, error_message = fetch_dna_from_genbank(accession_id)
        if sequence:
            st.success("Secuencia de ADN obtenida exitosamente!")
            visualize_dna(sequence)
        else:
            st.error(error_message)
    else:
        st.warning("Por favor, ingresa un Accession ID válido.")
