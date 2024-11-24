from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.SeqUtils import molecular_weight, gc_content
import streamlit as st
import matplotlib.pyplot as plt
from collections import Counter
from PIL import Image
import py3Dmol  # pip install py3Dmol==2.0.0.post2
import pandas as pd
import requests  # # python -m pip install requests
import io  # Import the 'io' module for StringIO
from io import StringIO
from io import BytesIO
import streamlit.components.v1 as components

# Creamos una barra
st.sidebar.header("Nucleótidos Operaciones 🧬")
sidebar_render = st.sidebar.radio("Opciones : ", ["Inicio", "Análisis de secuencia", "Parámetros de la estructura", "Secuencia de nucleótidos", "Visualizador de nucleótidos"])

# Página principal
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática: Análisis de Nucleótidos**')

    # Estilo de texto y colores
    st.markdown("""
    <style>
        .main-title {
            color: #4CAF50;  /* verde claro */
            font-size: 40px;
            font-weight: bold;
            text-align: center;
        }

        .text-block {
            color: #88dd9f;  /* verde más claro */
            font-size: 18px;
        }
        .team {
            font-style: bold;
            font-size: 16px;
            color: #7b9edd;  /* azul */
        }
    </style>
    """, unsafe_allow_html=True)

    # Título
    st.markdown('<div class="main-title">Bienvenido al Análisis de Nucleótidos</div>', unsafe_allow_html=True)

    # Descripción y subsecciones
    st.markdown("""
    <div class="text-block">
        Este tablero tiene el objetivo de facilitar el análisis y visualización de secuencias de nucleótidos como ADN y ARN.
        Explora diferentes herramientas interactivas para estudiar sus propiedades y estructura. Las secciones disponibles son:
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    - **🔬 Análisis de secuencia**: Carga archivos FASTA y analiza las secuencias de nucleótidos. Extrae información relevante como la composición de bases y propiedades biofísicas.
    - **🧬 Parámetros de la estructura**: Analiza características estructurales, como la cantidad de pares de bases, el contenido de GC y la estabilidad de la secuencia de ADN o ARN.
    - **🔍 Secuencia de nucleótidos**: Visualiza la secuencia y la proporción de bases de diversas secuencias de ADN/ARN.
    - **🌐 Visualización 3D de nucleótidos**: Introduce un código PDB y explora la estructura tridimensional de secuencias de ADN/ARN en modelos interactivos.
    """, unsafe_allow_html=True)

    # Línea divisoria
    st.markdown("<hr style='border:1px solid #ccc;'/>", unsafe_allow_html=True)

    # Mensaje motivador
    st.markdown("""
    <div class="text-block">
        ¡Explora las herramientas del lado izquierdo y haz un análisis profundo de las secuencias de nucleótidos que te interesen!
    </div>
    """, unsafe_allow_html=True)

    # Información del equipo
    st.markdown("<hr style='border:1px solid #ccc;'/>", unsafe_allow_html=True)
    st.markdown("<div class='team'>Equipo:</div>", unsafe_allow_html=True)
    st.markdown("""
    - **Camila García Rascón**
    - **Valeria Jara Salomón**
    """, unsafe_allow_html=True)

# Análisis de secuencia de nucleótidos
if sidebar_render == "Análisis de secuencia":
    st.title("🔬 Análisis de Secuencia de Nucleótidos")
    st.markdown("Sube tu secuencia de ADN o ARN en formato **FASTA** para analizarla ⬇️")

    # Función para leer y decodificar archivo FASTA
    def read_fasta_file(uploaded_file):
        fasta_file_content = uploaded_file.read().decode()
        return list(SeqIO.parse(StringIO(fasta_file_content), "fasta"))

    # Función para mostrar el contenido del archivo FASTA
    def display_fasta_file(sequence_contents):
        if len(sequence_contents) == 0:
            st.error("⚠️ Lo sentimos, no se encontraron secuencias en el archivo 😟")
        else:
            for i, record in enumerate(sequence_contents, start=1):
                with st.expander(f"🧬 Secuencia {i}: {record.id}"):
                    st.markdown(f"**🔖 ID:** `{record.id}`")
                    st.markdown(f"**🧾 Descripción:** {record.description}")
                    st.markdown(f"**🧪 Secuencia de nucleótidos:**")
                    st.code(str(record.seq), language="text")
                    st.markdown(f"**📏 Longitud de la secuencia:** `{len(record.seq)}`")
                st.divider()  # Línea divisoria entre secuencias

    # Subir archivo FASTA
    uploaded_file = st.file_uploader("📂 Sube tu archivo FASTA", type=["fasta"], help="Solo se admiten archivos con extensión .fasta")

    if uploaded_file is not None:
        with st.spinner("Procesando archivo... 🕒"):
            sequence_contents = read_fasta_file(uploaded_file)
            st.success("¡Archivo procesado exitosamente! 🎉", icon="✅")
        display_fasta_file(sequence_contents)
        st.info("Puedes copiar la secuencia para realizar otras operaciones.", icon="✨")


# Parámetros de la estructura de nucleótidos
if sidebar_render == "Parámetros de la estructura":
    st.title("🔬 Cálculos de parámetros de la estructura de Nucleótidos")
    st.markdown("Introduce tu secuencia de ADN o ARN y ajusta los valores necesarios para analizar sus propiedades estructurales. 🌟")

    # Entrada de secuencia
    sequence_input = st.text_area("✍️ Ingresa la secuencia de nucleótidos:")

    if st.button("⚡ ¡Calcular!"):
        if not sequence_input:
            st.error("Por favor, ingresa una secuencia para calcular sus propiedades.")
        else:
            # Análisis de la secuencia
            sequence_reference = Seq(str(sequence_input))
            st.markdown("<h3 style='text-align: center; color: #4CAF50;'>✨ Propiedades calculadas de la secuencia ✨</h3>", unsafe_allow_html=True)

            # Contenido de GC
            gc_percentage = round(gc_content(sequence_reference) * 100, 2)
            st.markdown("**1️⃣ Contenido de GC:**")
            st.info(f"⚖️ **Contenido de GC:** `{gc_percentage}%`")
            st.markdown("El contenido de GC es la proporción de bases de guanina y citosina en una secuencia de ADN.")

            # Peso molecular
            seq_weight = round(molecular_weight(sequence_reference), 2)
            st.markdown("**2️⃣ Peso molecular de la secuencia:**")
            st.info(f"⚖️ **Peso molecular:** `{seq_weight} Da`")
            st.markdown("El peso molecular de una secuencia de ADN o ARN es el valor calculado sumando el peso de todos los nucleótidos.")

            # Gráfico de la distribución de las bases
            base_counts = dict(Counter(sequence_reference))
            st.markdown("**3️⃣ Distribución de bases (A, T, C, G):**")
            st.bar_chart(base_counts)

            st.markdown("La distribución de bases muestra cuántas veces aparece cada base (A, T, C, G) en la secuencia.")

