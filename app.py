from Bio import Entrez, SeqIO
from Bio.SeqUtils import molecular_weight, gc_fraction
import streamlit as st
import matplotlib.pyplot as plt
from collections import Counter
import py3Dmol  # Para la visualización 3D
import pandas as pd

# Función para obtener la secuencia desde GenBank
def get_sequence_from_genbank(genbank_id):
    Entrez.email = "your_email@example.com"  # Cambia esto por tu correo
    try:
        # Buscar el ID en GenBank
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        return record
    except Exception as e:
        st.error(f"Error al obtener la secuencia de GenBank: {e}")
        return None

# Función para visualizar la secuencia de nucleótidos en 3D (py3Dmol)
def visualize_3D_sequence(sequence):
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(str(sequence), "fasta")
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    viewer.show()

# Página principal
st.sidebar.header("Nucleótidos Operaciones 🧬")
sidebar_render = st.sidebar.radio("Opciones : ", ["Inicio", "Análisis de secuencia de GenBank", "Parámetros estructurales de nucleótidos", "Visualización 3D de nucleótidos", "Distribución de bases", "Análisis de codones"])

# Página de inicio
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática: Análisis de Nucleótidos desde GenBank**')

    st.markdown("""
    Este tablero tiene el objetivo de facilitar el análisis y visualización de secuencias de nucleótidos como ADN y ARN. 
    Puedes ingresar un ID de GenBank para obtener la secuencia asociada y estudiar sus propiedades. Las secciones disponibles son:
    - **🔬 Análisis de secuencia de GenBank**: Ingresa un ID de GenBank y analiza la secuencia de nucleótidos.
    - **🧬 Parámetros estructurales de nucleótidos**: Analiza características como el contenido de GC y el peso molecular de las secuencias.
    - **🔍 Visualización 3D de nucleótidos**: Visualiza la secuencia de nucleótidos en 3D y realiza un análisis interactivo.
    - **📊 Distribución de bases**: Muestra un gráfico interactivo con la distribución de las bases de la secuencia.
    - **🔬 Análisis de codones**: Analiza los codones en la secuencia y su frecuencia.
    """)

# Análisis de secuencia de nucleótidos desde GenBank
if sidebar_render == "Análisis de secuencia de GenBank":
    st.title("🔬 Análisis de Secuencia de Nucleótidos desde GenBank")
    st.markdown("Ingresa un **ID de GenBank** para analizar la secuencia de ADN/ARN ⬇️")

    # Entrada de ID de GenBank
    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Mostrar el nombre del organismo
                organism = record.annotations.get("organism", "Desconocido")
                st.markdown(f"**🦠 Organismo o especie:** {organism}")

                # Mostrar la secuencia de nucleótidos
                st.markdown(f"**🔖 ID de GenBank:** `{record.id}`")
                st.markdown(f"**📜 Descripción:** {record.description}")
                st.markdown("**🧪 Secuencia de nucleótidos:**")
                st.code(str(record.seq), language="text")

                # Propiedades de la secuencia
                st.markdown("**1️⃣ Contenido de GC:**")
                gc_percentage = round(gc_fraction(record.seq) * 100, 2)
                st.info(f"⚖️ **Contenido de GC:** `{gc_percentage}%`")

                # Peso molecular de la secuencia
                seq_weight = round(molecular_weight(record.seq), 2)
                st.markdown("**2️⃣ Peso molecular de la secuencia:**")
                st.info(f"⚖️ **Peso molecular:** `{seq_weight} Da`")

                # Gráfico de la distribución de las bases
                base_counts = dict(Counter(record.seq))
                st.markdown("**3️⃣ Distribución de bases (A, T, C, G):**")
                st.bar_chart(base_counts)

# Parámetros estructurales de nucleótidos
if sidebar_render == "Parámetros estructurales de nucleótidos":
    st.title("🧬 Parámetros de la Estructura de Nucleótidos")
    st.markdown("Ingresa un **ID de GenBank** para calcular los parámetros de la secuencia de nucleótidos ⬇️")

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")
                
                # Mostrar nombre del organismo
                organism = record.annotations.get("organism", "Desconocido")
                st.markdown(f"**🦠 Organismo o especie:** {organism}")

                # Calcular y mostrar parámetros estructurales
                gc_percentage = round(gc_fraction(record.seq) * 100, 2)
                st.markdown(f"**1️⃣ Contenido de GC:** {gc_percentage}%")

                seq_weight = round(molecular_weight(record.seq), 2)
                st.markdown(f"**2️⃣ Peso molecular:** {seq_weight} Da")

                # Gráfico de la distribución de bases
                base_counts = dict(Counter(record.seq))
                st.markdown("**3️⃣ Distribución de bases (A, T, C, G):**")
                st.bar_chart(base_counts)

# Visualización 3D de nucleótidos
if sidebar_render == "Visualización 3D de nucleótidos":
    st.title("🌐 Visualización 3D de Nucleótidos")
    st.markdown("Ingresa un **ID de GenBank** para visualizar la secuencia de nucleótidos en 3D ⬇️")

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank para ver la secuencia 3D:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Visualizar la secuencia en 3D
                st.markdown(f"**🔖 ID de GenBank:** `{record.id}`")
                st.markdown(f"**📜 Descripción:** {record.description}")
                st.markdown("**🧪 Secuencia de nucleótidos:**")
                st.code(str(record.seq), language="text")

                # Visualizar la secuencia de nucleótidos en 3D
                visualize_3D_sequence(record.seq)

# Distribución de bases
if sidebar_render == "Distribución de bases":
    st.title("📊 Distribución de Bases")
    st.markdown("Muestra un gráfico interactivo de la distribución de bases (A, T, C, G) en la secuencia de nucleótidos ⬇️")

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Distribución de las bases (A, T, C, G)
                base_counts = dict(Counter(record.seq))
                st.markdown("**Distribución de bases (A, T, C, G):**")
                st.bar_chart(base_counts)
            
