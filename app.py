from Bio.Seq import Seq
from Bio import Entrez, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Blast import NCBIWWW, NCBIXML
import streamlit as st
from stmol import showmol
from stmol import *  # pip install stmol==0.0.9 , pip install ipython_genutils
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
st.sidebar.header("Proteínas Operaciones 🧬")
sidebar_render = st.sidebar.radio("Opciones : ", ["Inicio", "Análisis de secuencia", "Parámetros de la estructura", "Secuencia de aminoácidos de proteínas", "Visualizador de proteínas"])

# Página principal
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática: Análisis de Proteínas**')

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
    st.markdown('<div class="main-title">Bienvenido al Análisis de Proteínas</div>', unsafe_allow_html=True)

    # Descripción y subsecciones
    st.markdown("""
    <div class="text-block">
        Este tablero tiene el objetivo de facilitar el análisis y visualización de proteínas a partir de sus secuencias y estructuras. 
        Explora diferentes herramientas interactivas para estudiar sus propiedades y estructura. Las secciones disponibles son:
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    - **🔬 Análisis de secuencia**: Carga archivos FASTA y analiza las secuencias de proteínas. Extrae información relevante como la composición de aminoácidos y propiedades biofísicas.
    - **🧬 Parámetros de la estructura**: Calcula características estructurales, como el peso molecular, el punto isoeléctrico y la estabilidad de las proteínas, con un análisis detallado a nivel molecular.
    - **🔍 Secuencia de aminoácidos de proteínas**: Visualiza la secuencia y la proporción de átomos de diversas proteínas, con gráficos que permiten una mejor interpretación de sus características.
    - **🌐 Visualización 3D de proteínas**: Introduce un código PDB y explora la estructura tridimensional de proteínas en modelos interactivos. Personaliza la visualización y observa la estructura desde diferentes perspectivas.
    """, unsafe_allow_html=True)

    # Línea divisoria
    st.markdown("<hr style='border:1px solid #ccc;'/>", unsafe_allow_html=True)

    # Mensaje motivador
    st.markdown("""
    <div class="text-block">
        ¡Explora las herramientas del lado izquierdo y haz un análisis profundo de las proteínas que te interesen!
    </div>
    """, unsafe_allow_html=True)

    # Información del equipo
    st.markdown("<hr style='border:1px solid #ccc;'/>", unsafe_allow_html=True)
    st.markdown("<div class='team'>Equipo:</div>", unsafe_allow_html=True)
    st.markdown("""
    - **Camila García Rascón**
    - **Valeria Jara Salomón**
    """, unsafe_allow_html=True)

# Creamos Análisis de Secuencia
if sidebar_render == "Análisis de secuencia":
    st.title("🔬 Análisis de Secuencia")
    st.markdown("Sube tu secuencia de proteína en formato **FASTA** para analizarla ⬇️")

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
                    st.markdown(f"**🧪 Secuencia de aminoácidos:**")
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


# Parámetros de la estructura de proteínas
if sidebar_render == "Parámetros de la estructura":
    st.title("🔬 Cálculos de parámetros de la estructura de proteínas")
    st.markdown("Introduce tu secuencia de proteína y ajusta los valores necesarios para analizar sus propiedades estructurales. 🌟")

    # Entrada de secuencia y nivel de pH
    sequence_input = st.text_area("✍️ Ingresa la secuencia de aminoácidos:")
    pH = st.number_input("🌡️ ¿Con qué nivel de pH deseas analizar tu proteína?", min_value=0.0, max_value=14.0, value=7.0, step=0.1)
    
    if st.button("⚡ ¡Calcular!"):
        if not sequence_input:
            st.error("Por favor, ingresa una secuencia para calcular sus propiedades.")
        else:
            # Análisis de la proteína
            sequence_reference = ProteinAnalysis(str(sequence_input))
            st.markdown("<h3 style='text-align: center; color: #4CAF50;'>✨ Propiedades calculadas de la proteína ✨</h3>", unsafe_allow_html=True)

            # Número de aminoácidos
            st.markdown("**1️⃣ Número de aminoácidos:**")
            st.info(f"🔢 **Valor:** `{sequence_reference.count_amino_acids()}`")
            st.markdown("Los aminoácidos son moléculas que se combinan para formar proteínas. Los aminoácidos y las proteínas son los pilares fundamentales de la vida. Cuando las proteínas se digieren o se descomponen, el resultado son los aminoácidos.")

            # Peso molecular
            molecular_weight = round(sequence_reference.molecular_weight(), 2)
            st.markdown("**2️⃣ Peso molecular:**")
            st.info(f"⚖️ **Peso molecular:** `{molecular_weight} Da`")
            st.markdown("Los marcadores de peso molecular, o ladders, son un conjunto de estándares que se utilizan para determinar el tamaño aproximado de una proteína o un de fragmento de ácido nucleico procesado en un gel de electroforesis.")

            # Aromaticidad con barra de progreso
            aromaticity = round(sequence_reference.aromaticity(), 2)
            st.markdown("**3️⃣ Aromaticidad:**")
            st.info("Proporción de aminoácidos aromáticos en la proteína.")
            progress_value = int(aromaticity * 100)  # Convertir a porcentaje
            st.progress(progress_value)  # Barra de progreso de 0 a 100%
            st.markdown(f"**{progress_value}%** de la proteína tiene aminoácidos aromáticos.")
            st.markdown("La aromaticidad es una propiedad de las estructuras cíclicas, no saturadas, cuya estabilidad es superior a la de las estructuras de cadena abierta con igual número de enlaces múltiples.")

            # Índice de inestabilidad con barra visual
            instability_index = round(sequence_reference.instability_index(), 2)
            stability = "La proteína es inestable" if instability_index >= 40 else "La proteína es estable"
            st.markdown("**4️⃣ Índice de inestabilidad:**")
            st.info(f"📉 **Valor:** `{instability_index}` - ⚖️ **Estabilidad:** {stability}")
            progress_value_instability = min(int(instability_index), 100)  # Convertir a porcentaje y limitar a 100
            st.progress(min(int(instability_index), 100))  # Barra de
