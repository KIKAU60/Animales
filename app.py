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
    Entrez.email = "your_email@example.com"  # Cambia esto por tu correo
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

# Visualización 3D de ADN
if sidebar_render == "Visualización 3D de ADN":
    st.title("🌐 Visualización 3D de ADN")
    st.markdown("Ingresa un **ID de GenBank** para visualizar la secuencia de ADN en 3D ⬇️")

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank para ver la secuencia 3D:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Visualizar la secuencia en 3D
                st.markdown(f"**🔖 ID de GenBank:** `{record.id}`")
                st.markdown(f"**📜 Descripción:** {record.description}")
                st.markdown("**🧪 Secuencia de ADN:**")
                st.code(str(record.seq), language="text")

                # Visualizar la secuencia de ADN en 3D
                visualize_3D_dna(record.seq)

# Cantidad de proteínas codificadas, genes y cromosomas
if sidebar_render == "Cantidad de proteínas codificadas, genes y cromosomas":
    st.title("🧬 Cantidad de Proteínas Codificadas, Genes y Cromosomas")
    st.markdown("Ingresa un **ID de GenBank** para obtener la cantidad de proteínas codificadas, genes y cromosomas ⬇️")

    genbank_id = st.text_input("🧬 Ingresa el ID de GenBank:", "")

    if genbank_id:
        with st.spinner("Cargando datos desde GenBank... 🕒"):
            record = get_sequence_from_genbank(genbank_id)
            if record:
                st.success("¡Secuencia obtenida exitosamente! 🎉", icon="✅")

                # Obtener información adicional
                cds_count = len([feature for feature in record.features if feature.type == 'CDS'])
                gene_count = len([feature for feature in record.features if feature.type == 'gene'])
                chromosome_count = len([feature for feature in record.features if feature.type == 'chromosome'])

                # Mostrar los resultados
                st.markdown(f"**Proteínas Codificadas (CDS):** {cds_count}")
                st.markdown(f"**Cantidad de Genes:** {gene_count}")
                st.markdown(f"**Cantidad de Cromosomas:** {chromosome_count}")

                # Opcionalmente, podemos agregar imágenes ilustrativas si las tenemos disponibles
                st.image("proteins_coding.png", caption="Proteínas Codificadas", use_column_width=True)
                st.image("genes_count.png", caption="Cantidad de Genes", use_column_width=True)
                st.image("chromosomes_count.png", caption="Cantidad de Cromosomas", use_column_width=True

# Información adicional o conclusiones
if sidebar_render != "Inicio":
    st.sidebar.markdown("""
    Para obtener más detalles sobre cómo interpretar los resultados o cómo funciona el análisis de secuencias de GenBank, consulta la documentación de Biopython o el sitio web oficial de GenBank.
    Si deseas realizar otro análisis, simplemente elige una de las opciones en el menú lateral.
    """)

