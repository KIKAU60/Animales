from Bio import Entrez, SeqIO
from Bio.SeqUtils import molecular_weight, gc_fraction
import streamlit as st
import matplotlib.pyplot as plt
from collections import Counter
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

# Página principal
st.sidebar.header("Nucleótidos Operaciones 🧬")
sidebar_render = st.sidebar.radio("Opciones : ", ["Inicio", "Análisis de secuencia de GenBank", "Parámetros de la estructura", "Visualizador de nucleótidos"])

# Página de inicio
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática: Análisis de Nucleótidos desde GenBank**')

    # Descripción y subsecciones
    st.markdown("""
    Este tablero tiene el objetivo de facilitar el análisis y visualización de secuencias de nucleótidos como ADN y ARN. 
    Puedes ingresar un ID de GenBank para obtener la secuencia asociada y estudiar sus propiedades. Las secciones disponibles son:
    - **🔬 Análisis de secuencia de GenBank**: Ingresa un ID de GenBank y analiza la secuencia de nucleótidos.
    - **🧬 Parámetros de la estructura**: Analiza características estructurales como el contenido de GC y el peso molecular de las secuencias.
    - **🔍 Secuencia de nucleótidos**: Visualiza la secuencia y la proporción de bases de diversas secuencias de ADN/ARN.
    - **🌐 Visualización 3D de nucleótidos**: Introduce un código PDB y explora la estructura tridimensional de secuencias de ADN/ARN.
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
                gc_percentage = round(gc_fraction(record.seq) * 100, 2)  # Usar gc_fraction en lugar de gc_content
                st.info(f"⚖️ **Contenido de GC:** `{gc_percentage}%`")
                st.markdown("El contenido de GC es la proporción de bases de guanina y citosina en una secuencia de ADN.")
                
                # Peso molecular de la secuencia
                seq_weight = round(molecular_weight(record.seq), 2)
                st.markdown("**2️⃣ Peso molecular de la secuencia:**")
                st.info(f"⚖️ **Peso molecular:** `{seq_weight} Da`")
                st.markdown("El peso molecular de una secuencia de ADN o ARN es el valor calculado sumando el peso de todos los nucleótidos.")

                # Gráfico de la distribución de las bases
                base_counts = dict(Counter(record.seq))
                st.markdown("**3️⃣ Distribución de bases (A, T, C, G):**")
                st.bar_chart(base_counts)

                st.markdown("La distribución de bases muestra cuántas veces aparece cada base (A, T, C, G) en la secuencia.")

# Parámetros de la estructura de nucleótidos desde GenBank
if sidebar_render == "Parámetros de la estructura":
    st.title("🧬 Parámetros de la Estructura de Nucleótidos desde GenBank")
    st.markdown("Ingresa un **ID de GenBank** para calcular los parámetros de la secuencia de nucleótidos ⬇️")

    # Entrada de ID de GenBank
    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank para calcular sus parámetros:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Mostrar el nombre del organismo
                organism = record.annotations.get("organism", "Desconocido")
                st.markdown(f"**🦠 Organismo o especie:** {organism}")

                # Contenido de GC
                gc_percentage = round(gc_fraction(record.seq) * 100, 2)
                st.markdown("**1️⃣ Contenido de GC:**")
                st.info(f"⚖️ **Contenido de GC:** `{gc_percentage}%`")

                # Peso molecular
                seq_weight = round(molecular_weight(record.seq), 2)
                st.markdown("**2️⃣ Peso molecular:**")
                st.info(f"⚖️ **Peso molecular:** `{seq_weight} Da`")

                # Distribución de las bases
                base_counts = dict(Counter(record.seq))
                st.markdown("**3️⃣ Distribución de bases (A, T, C, G):**")
                st.bar_chart(base_counts)

                st.markdown("La distribución de bases muestra cuántas veces aparece cada base (A, T, C, G) en la secuencia.")

# Visualizador de nucleótidos desde GenBank
if sidebar_render == "Visualizador de nucleótidos":
    st.title("🌐 Visualizador de Nucleótidos desde GenBank")
    st.markdown("Ingresa un **ID de GenBank** para visualizar la secuencia de nucleótidos ⬇️")

    # Entrada de ID de GenBank
    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank para ver la secuencia:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Mostrar la secuencia
                st.markdown(f"**🔖 ID de GenBank:** `{record.id}`")
                st.markdown(f"**📜 Descripción:** {record.description}")
                st.markdown("**🧪 Secuencia de nucleótidos:**")
                st.code(str(record.seq), language="text")

                # Visualizar la distribución de bases
                base_counts = dict(Counter(record.seq))
                st.markdown("**🔢 Distribución de bases (A, T, C, G):**")
                st.bar_chart(base_counts)

                st.markdown("Esta es la distribución de bases de la secuencia de ADN o ARN.")
