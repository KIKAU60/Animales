import streamlit as st
from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt

# Configurar Entrez
Entrez.email = "a223201128@unison.mx"  # Asegúrate de poner tu correo electrónico

# Función para obtener la secuencia de ADN y datos del gen desde GenBank
def fetch_gene_data_from_genbank(accession_id):
    try:
        # Buscar el archivo GenBank utilizando el ID de acceso
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Verificar si la secuencia es válida
        if hasattr(record, 'seq') and record.seq:
            # Extraer información adicional del GenBank (como características)
            gene_description = record.description
            gene_features = record.features
            organism = record.annotations.get("organism", "Desconocido")
            sequence = record.seq
            return sequence, gene_description, gene_features, organism
        else:
            return None, "Secuencia no encontrada en el archivo de GenBank.", None, None

    except Exception as e:
        return None, f"Error al obtener la secuencia de GenBank: {e}", None, None

# Función para calcular el contenido GC de la secuencia
def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in 'GCgc')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

# Función para mostrar la información del gen
def visualize_gene_info(sequence, gene_description, gene_features, organism):
    if sequence and len(sequence) > 0:
        sequence = sequence.upper()  # Convertir la secuencia a mayúsculas
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)

        # Mostrar la descripción del gen y organismo
        st.subheader("Descripción del Gen:")
        st.write(gene_description)
        st.subheader("Organismo:")
        st.write(organism)

        # Mostrar información adicional
        st.subheader("Información sobre el Gen:")
        st.write(f"Longitud de la secuencia: {length} bases")
        st.write(f"Contenido GC: {gc_content:.2f}%")
        
        # Mostrar características del gen
        st.subheader("Características del Gen:")
        for feature in gene_features:
            if feature.type == "CDS":  # CDS (Código de secuencia) es una característica común de genes
                st.write(f"- {feature.type} en las posiciones {feature.location}")
                if "product" in feature.qualifiers:
                    st.write(f"  Función: {feature.qualifiers['product'][0]}")

        # Mostrar gráfico de contenido GC
        st.subheader("Gráfico de distribución de bases")
        fig, ax = plt.subplots()
        ax.bar(["A", "T", "G", "C"], [sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C')])
        ax.set_ylabel("Número de bases")
        ax.set_title("Distribución de bases en la secuencia de ADN")
        st.pyplot(fig)

    else:
        st.warning("La secuencia de ADN no se pudo obtener o está vacía. Por favor, verifica el ID de acceso.")

# Título de la aplicación
st.title("Análisis de Gen de Animal desde GenBank")

# Entrada para el ID de GenBank
accession_id = st.text_input("Introduce el ID de acceso de GenBank (Accession ID):", 
                             help="Puedes buscar el ID de acceso de un gen o genoma en GenBank (por ejemplo, NC_000913)")

# Botón para obtener la secuencia desde GenBank
if st.button("Obtener información genética desde GenBank"):
    if accession_id:
        sequence, gene_description, gene_features, organism = fetch_gene_data_from_genbank(accession_id)
        if sequence:
            st.success("Secuencia de ADN obtenida exitosamente!")
            visualize_gene_info(sequence, gene_description, gene_features, organism)
        else:
            st.error(gene_description)
    else:
        st.warning("Por favor, ingresa un Accession ID válido.")
