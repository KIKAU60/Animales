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
    st.title('🧬 **Bioinformática: Análisis de Nucleótidos desde GenBank**')

    st.markdown(""" 
    Este tablero tiene el objetivo de facilitar el análisis y visualización de secuencias de nucleótidos como ADN y ARN. 
    Puedes ingresar un ID de GenBank para obtener la secuencia asociada y estudiar sus propiedades. Las secciones disponibles son:
    - **🔬 Frecuencia de codones**: Análisis de la frecuencia de codones en la secuencia de nucleótidos.
    - **📊 Distribución de bases nitrogenadas**: Analiza la distribución de las bases nitrogenadas A, T, C, G de la secuencia.
    - **🔬 Análisis de Motivos Conservados**: Análisis de secuencias de ADN conservadas en todo el genoma.
    - **🔬 Cálculo de Enriquecimiento de GC**: Calcula y visualiza el contenido de GC a lo largo de la secuencia.
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
    st.markdown("Ingresa un **ID de GenBank** para analizar la distribución de bases A, T, C, G ⬇️")

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
    st.markdown("Introduce el ID de GenBank para analizar los motivos conservados en la secuencia de ADN. 🌟")

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
    st.markdown("Introduce el ID de GenBank para analizar el contenido de GC en la secuencia de ADN. 🌟")

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
