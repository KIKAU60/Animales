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

# Función para visualizar la secuencia de ADN en 3D (py3Dmol)
def visualize_3D_dna(sequence):
    viewer = py3Dmol.view(width=800, height=600)
    viewer.addModel(str(sequence), "fasta")
    viewer.setStyle({'stick': {}})
    viewer.setBackgroundColor('white')
    viewer.zoomTo()
    viewer.show()

# Página principal
st.sidebar.header("Nucleótidos Operaciones 🧬")
sidebar_render = st.sidebar.radio("Opciones : ", ["Inicio", "Frecuencia de codones", "Distribución de bases nitrogenadas", "Visualización 3D de ADN", "Cantidad de proteínas codificadas, genes y cromosomas"])

# Página de inicio
if sidebar_render == "Inicio":
    st.title('🧬 **Bioinformática: Análisis de Nucleótidos desde GenBank**')

    st.markdown("""
    Este tablero tiene el objetivo de facilitar el análisis y visualización de secuencias de nucleótidos como ADN y ARN. 
    Puedes ingresar un ID de GenBank para obtener la secuencia asociada y estudiar sus propiedades. Las secciones disponibles son:
    - **🔬 Frecuencia de codones**: Análisis de la frecuencia de codones en la secuencia de nucleótidos.
    - **📊 Distribución de bases nitrogenadas**: Analiza la distribución de las bases nitrogenadas A, T, C, G de la secuencia.
    - **🌐 Visualización 3D de ADN**: Visualiza la secuencia de ADN en 3D.
    - **🧬 Cantidad de proteínas codificadas, genes y cromosomas**: Muestra la cantidad de proteínas codificadas, genes y cromosomas, ilustrados.
    """)

# Frecuencia de codones
if sidebar_render == "Frecuencia de codones":
    st.title("🔬 Frecuencia de Codones en la Secuencia de Nucleótidos")
    st.markdown("Ingresa un **ID de GenBank** para analizar la frecuencia de codones ⬇️")

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

def fetch_genbank_record(genbank_id):
    """
    Esta función obtiene el registro de GenBank usando el ID proporcionado.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        return record
    except Exception as e:
        st.error(f"Error al recuperar el ID de GenBank: {e}")
        return None

def fetch_genbank_record(genbank_id):
    """
    Esta función obtiene el registro de GenBank usando el ID proporcionado.
    """
    try:
        handle = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        return record
    except Exception as e:
        st.error(f"Error al recuperar el ID de GenBank: {e}")
        return None

# 1. Análisis de Motivos Regulatorios en el ADN
if sidebar_render == "Análisis de Motivos Regulatorios":
    st.title("🔬 Análisis de Motivos Regulatorios en el ADN")
    st.markdown("Introduce el ID de GenBank para analizar los motivos regulatorios en el ADN. 🌟")

    # Entrada para el ID de GenBank
    genbank_id = st.text_input("✍️ Ingresa el ID de GenBank", "NM_001301717")  # ID de ejemplo

    if st.button("⚡ ¡Analizar!"):
        if not genbank_id:
            st.error("Por favor, ingresa un ID de GenBank válido.")
        else:
            with st.spinner("Cargando información desde GenBank... 🕒"):
                # Acceder al registro GenBank con Biopython
                record = fetch_genbank_record(genbank_id)
                if record:
                    # Obtener la secuencia de ADN
                    sequence = record.seq
                    
                    # Simulación de análisis de motivos regulatorios (esto es un ejemplo)
                    # Supón que hemos identificado algunos motivos regulatorios
                    regulatory_motifs = ['TATA', 'CAAT', 'GC-box']
                    motif_positions = [i for i in range(len(sequence)) if sequence[i:i+4] in regulatory_motifs]

                    # Visualización de los motivos regulatorios en un gráfico de barras
                    st.markdown("**🔬 Posiciones de Motivos Regulatorios**")
                    fig = go.Figure(data=[go.Bar(
                        x=list(range(len(motif_positions))),
                        y=[1]*len(motif_positions),  # Solo para ilustrar la presencia de los motivos
                        marker=dict(color='royalblue')
                    )])
                    fig.update_layout(
                        title="Posiciones de Motivos Regulatorios en la Secuencia",
                        xaxis_title="Posición en la secuencia",
                        yaxis_title="Presencia de Motivo",
                        template="plotly_dark"
                    )
                    st.plotly_chart(fig)

# 2. Análisis de SNPs (Polimorfismos de Nucleótido Único)
if sidebar_render == "Análisis de SNPs":
    st.title("🔬 Análisis de SNPs (Polimorfismos de Nucleótido Único)")
    st.markdown("Introduce el ID de GenBank para analizar los SNPs en la secuencia de ADN. 🌟")

    # Entrada para el ID de GenBank
    genbank_id = st.text_input("✍️ Ingresa el ID de GenBank", "NM_001301717")  # ID de ejemplo

    if st.button("⚡ ¡Analizar!"):
        if not genbank_id:
            st.error("Por favor, ingresa un ID de GenBank válido.")
        else:
            with st.spinner("Cargando información desde GenBank... 🕒"):
                # Acceder al registro GenBank con Biopython
                record = fetch_genbank_record(genbank_id)
                if record:
                    # Obtener la secuencia de ADN
                    sequence = record.seq
                    
                    # Simulación de análisis de SNPs (esto es un ejemplo)
                    # Supón que hemos identificado algunos SNPs
                    snps = [i for i in range(1, len(sequence)) if sequence[i] != sequence[i-1]]

                    # Visualización de los SNPs en un gráfico interactivo
                    st.markdown("**🔬 Posiciones de SNPs**")
                    fig = go.Figure(data=[go.Scatter(
                        x=snps,
                        y=[1]*len(snps),  # Solo para ilustrar las posiciones de los SNPs
                        mode='markers',
                        marker=dict(color='red', size=10)
                    )])
                    fig.update_layout(
                        title="Posiciones de SNPs en la Secuencia",
                        xaxis_title="Posición en la secuencia",
                        yaxis_title="SNP Detectado",
                        template="plotly_dark"
                    )
                    st.plotly_chart(fig)
