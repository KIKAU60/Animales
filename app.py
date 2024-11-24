import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
import re
from Bio import Entrez
from Bio import SeqIO

# Configuración de Entrez: Debes proporcionar tu correo para cumplir con los requisitos de GenBank
Entrez.email = "a223201128@unison.mx"  # Sustituye con tu correo electrónico

# Función para obtener la secuencia de ADN desde GenBank
def obtener_secuencia_genbank(accession_number):
    try:
        # Fetch la secuencia de GenBank usando el número de acceso
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")  # Lee el archivo de GenBank
        secuencia_adn = str(record.seq)  # Convierte la secuencia en cadena
        nombre_organismo = record.annotations.get('organism', 'Organismo desconocido')  # Obtener el nombre del organismo
        return secuencia_adn, nombre_organismo
    except Exception as e:
        st.error(f"Error al obtener la secuencia de GenBank: {e}")
        return None, None

# Función para calcular las proporciones de nucleótidos (A, T, C, G)
def calcular_proporcion_nucleotidos(secuencia):
    counter = Counter(secuencia)
    total = len(secuencia)
    proporciones = [counter['A'] / total, counter['T'] / total, counter['C'] / total, counter['G'] / total]
    return proporciones

# Función para graficar la frecuencia de codones (tripletas) usando Plotly
def graficar_codones_interactivo(secuencia):
    codones = [secuencia[i:i+3] for i in range(0, len(secuencia), 3)]
    counter = Counter(codones)
    codones_unicos = list(counter.keys())
    frecuencias = list(counter.values())

    # Gráfico interactivo usando Plotly
    fig = go.Figure([go.Bar(x=codones_unicos, y=frecuencias, marker_color='skyblue')])
    fig.update_layout(
        title="Frecuencia de Codones en la Secuencia de ADN",
        xaxis_title="Codón",
        yaxis_title="Frecuencia",
        xaxis_tickangle=-45
    )
    st.plotly_chart(fig)

# Función para analizar secuencias de replicación (orígenes de replicación)
def analizar_replicacion(secuencia):
    # Detectamos secuencias comunes en orígenes de replicación, como la secuencia "ATG" o similares
    secuencias_replicacion = re.findall(r'ATG[AGCT]{5,10}ATG', secuencia)  # Secuencias con características de replicación
    return secuencias_replicacion

# Función para generar la hélice de ADN en 3D usando Plotly
def mostrar_helice_3d():
    # Ejemplo simple de la hélice de ADN utilizando coordenadas 3D
    phi = np.linspace(0, 2 * np.pi, 100)
    z = np.linspace(-5, 5, 100)
    x = np.sin(phi) * np.exp(0.1 * z)
    y = np.cos(phi) * np.exp(0.1 * z)

    # Crear la figura 3D con Plotly
    fig = go.Figure(data=[go.Scatter3d(
        x=x, y=y, z=z,
        mode='lines',
        line=dict(color='blue', width=4)
    )])
    
    fig.update_layout(
        title="Estructura 3D de la Hélice de ADN",
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z'
        )
    )
    
    st.plotly_chart(fig)

# Función principal para la aplicación Streamlit
def main():
    # Barra lateral de navegación
    st.sidebar.title("Análisis de ADN")
    st.sidebar.markdown("""
    Bienvenido al análisis interactivo de secuencias de ADN.
    Elige una opción para ver el análisis correspondiente:
    """)
    
    # Selección en barra lateral
    opcion = st.sidebar.radio(
        "Selecciona una opción",
        ["Introducción", "Análisis de Secuencia", "Visualización Interactiva"]
    )

    if opcion == "Introducción":
        st.title("Bioinformática: Características del ADN animal")
        st.write("""
            Esta aplicación te permitirá analizar secuencias de ADN obtenidas desde GenBank, 
            y visualizar diferentes características e ilustraciones interactivas sobre las secuencias.
        """)
    
    elif opcion == "Análisis de Secuencia":
        st.title("Introduce el Número de Acceso de GenBank")
        accession_number = st.text_input("Introduce el número de acceso de GenBank (Ej. NM_001003222):")
        
        if accession_number:
            secuencia_adn, nombre_organismo = obtener_secuencia_genbank(accession_number)
            
            if secuencia_adn:
                st.write(f"Secuencia de ADN obtenida de GenBank (Accession Number: {accession_number}):")
                st.text(secuencia_adn)
                st.write(f"**Organismo:** {nombre_organismo}")
            else:
                st.write("No se pudo obtener la secuencia. Verifica el número de acceso.")
    
    elif opcion == "Visualización Interactiva":
        st.title("Selecciona una visualización")
        
        ilustracion = st.selectbox(
            "Selecciona la ilustración para la secuencia de ADN:",
            ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Secuencias de Replicación', 'Hélice 3D de ADN']
        )
        
        # Descripción e ilustración según la opción seleccionada
        if ilustracion == 'Proporciones de Nucleótidos':
            st.markdown("### Proporciones de Nucleótidos")
            st.write("""
                Esta visualización muestra la distribución de los cuatro nucleótidos (A, T, C, G) 
                en la secuencia de ADN introducida.
            """)
            if accession_number:
                proporciones = calcular_proporcion_nucleotidos(secuencia_adn)
                fig, ax = plt.subplots(figsize=(8, 8))
                ax.pie(proporciones, labels=['A', 'T', 'C', 'G'], autopct='%1.1f%%', startangle=140)
                ax.set_title("Proporción de Nucleótidos en la Secuencia de ADN")
                st.pyplot(fig)

        elif ilustracion == 'Frecuencia de Codones':
            st.markdown("### Frecuencia de Codones")
            st.write("""
                Esta visualización muestra la frecuencia de aparición de los codones (tripletas de nucleótidos) 
                en la secuencia de ADN.
            """)
            if accession_number:
                graficar_codones_interactivo(secuencia_adn)

        elif ilustracion == 'Secuencias de Replicación':
            st.markdown("### Secuencias de Replicación")
            st.write("""
                Aquí se detectan posibles secuencias que puedan indicar orígenes de replicación en la secuencia de ADN.
            """)
            if accession_number:
                secuencias_replicacion = analizar_replicacion(secuencia_adn)
                if secuencias_replicacion:
                    st.write("Se han encontrado las siguientes secuencias de replicación (orígenes de replicación):")
                    st.write(secuencias_replicacion)
                else:
                    st.write("No se encontraron secuencias de replicación en la secuencia.")

        elif ilustracion == 'Hélice 3D de ADN':
            st.markdown("### Hélice 3D de ADN")
            st.write("""
                Esta visualización genera una representación 3D de la hélice de ADN utilizando coordenadas 3D.
            """)
            mostrar_helice_3d()

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
