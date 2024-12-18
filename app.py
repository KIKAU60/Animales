from Bio import Entrez, SeqIO
from Bio.SeqUtils import molecular_weight, gc_fraction
import streamlit as st
import matplotlib.pyplot as plt
from collections import Counter
import py3Dmol  # Para la visualización 3D
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

# Función para obtener la secuencia desde GenBank
def get_sequence_from_genbank(genbank_id):
    Entrez.email = "a223201128@unison.mx"  # Cambia esto por tu correo
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
sidebar_render = st.sidebar.radio("Opciones : ", ["Inicio", "Frecuencia de codones", "Distribución de bases nitrogenadas", "Análisis de Motivos Conservados", "Cálculo de Enriquecimiento de GC"])

# Página de inicio
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática**')
    st.markdown("""
        ## Bienvenidos al tablero de análisis de secuencias nucleotídicas 🧬

        Este tablero está diseñado para facilitar el análisis y la visualización de secuencias de ADN y ARN. 
        A través de este portal, puedes ingresar un **ID de GenBank** para obtener la secuencia asociada y estudiar sus propiedades fundamentales.

        ### Secciones disponibles:

        - **🔬 Frecuencia de codones**: Analiza la frecuencia de codones presentes en la secuencia de nucleótidos.
        - **📊 Distribución de bases nitrogenadas**: Examina la distribución de las bases nitrogenadas (A, T, C, G) a lo largo de la secuencia.
        - **🔬 Análisis de Motivos Conservados**: Estudia los motivos conservados dentro de las secuencias de ADN a lo largo del genoma.
        - **🔬 Cálculo de Enriquecimiento de GC**: Calcula y visualiza el contenido de GC a lo largo de la secuencia.

        ### **Equipo**:
        - Ana Camila Gracia Barroso
        - Ana Paola Terán Rascón
        - Diana Lizeth Villaescusa Guillén

        ¡Disfruta del análisis y la exploración de las secuencias!
    """)

# Frecuencia de codones
if sidebar_render == "Frecuencia de codones":
    st.title("🔬 Frecuencia de Codones en la Secuencia de Nucleótidos")
    st.markdown("""
    **Introducción:**
    Los **codones** son secuencias de tres nucleótidos consecutivos en el ADN que especifican la incorporación de un aminoácido particular durante la traducción de proteínas. La **frecuencia de codones** se refiere a cuántas veces aparecen cada uno de los 64 codones posibles en una secuencia de ADN. Este análisis es crucial para comprender cómo una célula utiliza su código genético, ya que ciertos codones pueden ser más frecuentes que otros debido a su eficiencia en la traducción o su influencia en la estabilidad del ARNm.
    
    **Aplicaciones:**
    - Optimización de la expresión genética.
    - Estudio de genes en organismos.
    - Detección de mutaciones.
    """)

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank:", "")

    # Sugerir ejemplos si no se ingresa un ID o si se ingresa uno incorrecto
    if not genbank_id:
        st.markdown("""
            **Ejemplos de IDs de GenBank**:
            - `NM_001301717`: Gen humano relacionado con la distrofia muscular de Duchenne.
            - `NM_001301806`: Gen humano relacionado con el sistema inmune.
            - `X16067`: ADN mitocondrial de Homo sapiens.
            - `EU551120`: ADN ribosómico de Escherichia coli.
            - `AY597007`: ADN ribosómico 16S de Pseudomonas aeruginosa.
        """)

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Codificar la secuencia y obtener frecuencia de codones
                codons = [str(record.seq[i:i+3]) for i in range(0, len(record.seq), 3)]
                codon_counts = dict(Counter(codons))

                # Crear la gráfica de barras interactiva para la frecuencia de codones
                st.markdown("**Frecuencia de Codones:**")
                codon_df = pd.DataFrame(list(codon_counts.items()), columns=['Codón', 'Frecuencia'])
                fig = px.bar(codon_df, x='Codón', y='Frecuencia', title="Frecuencia de Codones en la Secuencia", labels={'Codón': 'Codón', 'Frecuencia': 'Frecuencia'})
                st.plotly_chart(fig)

# Distribución de bases nitrogenadas
if sidebar_render == "Distribución de bases nitrogenadas":
    st.title("📊 Distribución de Bases Nitrogenadas")
    st.markdown("""
    **Introducción:**
    El ADN está compuesto por cuatro tipos de bases nitrogenadas: **adenina (A)**, **timina (T)**, **citoína (C)** y **guanina (G)**. La **distribución de bases nitrogenadas** describe cómo se distribuyen estas bases a lo largo de la secuencia de ADN. Cada organismo tiene una distribución particular, que puede ser informativa sobre su estructura genética, evolución, y función biológica.
    
    **Aplicaciones:**
    - Identificación de regiones codificantes y no codificantes.
    - Estudio de genomas.
    - Análisis de regiones ricas en GC.
    """)

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Distribución de bases
                base_counts = dict(Counter(record.seq))

                # Crear gráfico de pastel interactivo para la distribución de bases nitrogenadas
                st.markdown("**Distribución de Bases Nitrogenadas (A, T, C, G):**")
                base_df = pd.DataFrame(list(base_counts.items()), columns=['Base', 'Cantidad'])
                fig = px.pie(base_df, names='Base', values='Cantidad', title="Distribución de Bases Nitrogenadas")
                fig.update_traces(textinfo="percent+label", pull=[0.1, 0.1, 0.1, 0.1])
                st.plotly_chart(fig)

# 1. Análisis de Motivos Conservados con Heatmap
if sidebar_render == "Análisis de Motivos Conservados":
    st.title("🔬 Análisis de Motivos Conservados")
    st.markdown("""
    **Introducción:**
    Los **motivos conservados** son secuencias de nucleótidos o aminoácidos que se mantienen sin cambios a lo largo de la evolución debido a su función biológica esencial. En el análisis de **motivos conservados**, buscamos identificar estos patrones repetitivos en secuencias de ADN, ya que pueden ser cruciales para funciones específicas como la unión de proteínas, la replicación del ADN o la transcripción.
    
    **Aplicaciones:**
    - Predicción de sitios funcionales.
    - Estudio de evolución molecular.
    - Desarrollo de fármacos.
    """)

    # Entrada para el ID de GenBank
    genbank_id = st.text_input("✍️ Ingresa el ID de GenBank", "NM_001301717")  # ID de ejemplo

    if st.button("⚡ ¡Analizar!"):
        if not genbank_id:
            st.error("Por favor, ingresa un ID de GenBank válido.")
        else:
            with st.spinner("Cargando información desde GenBank... 🕒"):
                # Acceder al registro GenBank con Biopython
                record = get_sequence_from_genbank(genbank_id)
                if record:
                    # Obtener la secuencia de ADN
                    sequence = record.seq
                    
                    # Supongamos que hemos identificado algunos motivos conservados
                    conserved_motifs = ['ATG', 'TAA', 'GGT']  # Motivos conservados de ejemplo
                    motif_positions = [i for i in range(len(sequence)) if sequence[i:i+3] in conserved_motifs]

                    # Convertimos las posiciones de los motivos conservados en una matriz de intensidad
                    motif_matrix = [0] * len(sequence)
                    for pos in motif_positions:
                        motif_matrix[pos] = 1  # Marca con 1 donde el motivo está presente

                    # Crear una visualización de Heatmap con la intensidad de los motivos conservados
                    st.markdown("**🔬 Heatmap de Motivos Conservados a lo largo de la Secuencia**")

                    # Crear el heatmap usando Plotly
                    fig = go.Figure(data=go.Heatmap(
                        z=[motif_matrix],  # Los datos de la matriz
                        colorscale='Viridis',  # Escala de colores (puedes elegir otras como 'Cividis', 'Inferno', etc.)
                        colorbar=dict(title="Intensidad"),
                        showscale=True,  # Muestra la escala de colores
                        zmin=0,  # Mínimo valor de la escala
                        zmax=1,  # Máximo valor de la escala
                    ))

                    # Mejora de la presentación del gráfico
                    fig.update_layout(
                        title="Heatmap de Motivos Conservados en la Secuencia de ADN",
                        xaxis_title="Posiciones en la Secuencia de ADN",
                        yaxis_title="Motivos Conservados (0 = Ausente, 1 = Presente)",
                        template="plotly_dark",  # Estilo oscuro
                        plot_bgcolor="black",
                        paper_bgcolor="rgb(17, 17, 17)",
                        showlegend=False,
                        xaxis=dict(showgrid=True, zeroline=False),
                        yaxis=dict(showgrid=True, zeroline=False),
                    )

                    # Mostrar el gráfico interactivo
                    st.plotly_chart(fig)

# Cálculo de Enriquecimiento de GC con gráfico de líneas
if sidebar_render == "Cálculo de Enriquecimiento de GC":
    st.title("🔬 Cálculo de Enriquecimiento de GC")
    st.markdown("""
    **Introducción:**
    El contenido de **GC** se refiere a la proporción de bases **guanina (G)** y **citoína (C)** en una secuencia de ADN. El **enriquecimiento de GC** describe cómo varía este contenido a lo largo de una secuencia, y puede tener implicaciones sobre la estabilidad y la estructura del ADN. Las secuencias ricas en GC suelen ser más estables debido a los enlaces triples que unen G y C, en comparación con los enlaces dobles entre A y T.
    
    **Aplicaciones:**
    - Estabilidad estructural.
    - Análisis de genomas.
    - Estudio de adaptación evolutiva.
    """)

    genbank_id = st.text_input("✍️ Ingresa el ID de GenBank", "NM_001301717")  # ID de ejemplo

    if st.button("⚡ ¡Analizar!"):
        if not genbank_id:
            st.error("Por favor, ingresa un ID de GenBank válido.")
        else:
            with st.spinner("Cargando información desde GenBank... 🕒"):
                # Acceder al registro GenBank con Biopython
                record = get_sequence_from_genbank(genbank_id)
                if record:
                    # Obtener la secuencia de ADN
                    sequence = record.seq

                    # Calcular el contenido de GC por bloques de 100 nucleótidos
                    gc_blocks = [gc_fraction(sequence[i:i+100]) * 100 for i in range(0, len(sequence), 100)]
                    block_indices = [i for i in range(len(gc_blocks))]

                    # Crear gráfico de líneas interactivo para el contenido de GC en bloques
                    st.markdown("**🔬 Enriquecimiento de GC en la Secuencia**")
                    fig = go.Figure(data=[go.Scatter(
                        x=block_indices,
                        y=gc_blocks,
                        mode='lines+markers',  # Línea con marcadores
                        marker=dict(color='green'),
                        line=dict(width=2)
                    )])

                    # Actualización del gráfico para mejorar la visualización
                    fig.update_layout(
                        title="Enriquecimiento de GC a lo largo de la secuencia",
                        xaxis_title="Bloque de nucleótidos",
                        yaxis_title="Contenido de GC (%)",
                        template="plotly_dark",
                        xaxis=dict(
                            tickmode='array',
                            tickvals=block_indices,
                            ticktext=[f"Bloque {i+1}" for i in block_indices],
                        ),
                        yaxis=dict(
                            range=[0, 100],
                            ticksuffix="%"
                        ),
                        plot_bgcolor="black",
                        paper_bgcolor="rgb(17, 17, 17)",
                    )

                    # Mostrar el gráfico de líneas interactivo
                    st.plotly_chart(fig)
