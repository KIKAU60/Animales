import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
import seaborn as sns
import random
import re

# Secuencias de ADN de animales
secuencias_adn = {
    "Elefante": "AGCTGACGTAGCGTACGTAAGCTGACTGA",
    "Perro": "ATCGAGCTGGTAGCGGATCGAAGTCTAGG",
    "Gato": "AAGGCTAGCTAGGTACGTCGAAGTCGAGT",
    "Caballo": "AGGTCGACGTTGAGTCTGAGTGAGTCGA",
    "León": "ATGCGATCGTACGAGTGTAGCTAGCGTA",
    "Tigre": "CTGAGTGAGTCGATAGCGATGCAGTCAG",
    "Delfín": "AGTCTGATCGGAGTCTACGAGAGTCTGA",
    "Ballena": "ACGTGAGTACGAGTGTACGTAGTGACTG",
    "Cebra": "ATGAGTCTAGGATCGAGTACGAGGTCTGA",
    "Rinoceronte": "AGTCGTAGGCTAGCTGACGTAGCGTAGC"
}

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

# Función para detectar secuencias de regulación
def detectar_secuencias_reguladoras(secuencia):
    # Supongamos que los promotores y elementos potenciadores tienen secuencias específicas
    secuencias_reguladoras = re.findall(r'TATAAA|CAAT|GC-box', secuencia)
    return secuencias_reguladoras

# Función principal para la aplicación Streamlit
def main():
    st.title("Análisis de Secuencias de ADN de Animales")
    
    # Selección de un animal
    animal = st.selectbox("Selecciona un animal:", list(secuencias_adn.keys()))
    
    # Obtener la secuencia de ADN del animal seleccionado
    secuencia_adn = secuencias_adn[animal]
    
    # Mostrar la secuencia seleccionada
    st.write(f"Secuencia de ADN del {animal}:")
    st.text(secuencia_adn)
    
    # Opciones de ilustración
    ilustracion = st.selectbox(
        "Selecciona la ilustración para la secuencia de ADN:",
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Secuencias de Replicación', 'Secuencias de Regulación']
    )
    
    # Actualizar visualización según la opción seleccionada
    if ilustracion == 'Proporciones de Nucleótidos':
        proporciones = calcular_proporcion_nucleotidos(secuencia_adn)
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.pie(proporciones, labels=['A', 'T', 'C', 'G'], autopct='%1.1f%%', startangle=140)
        ax.set_title("Proporción de Nucleótidos en la Secuencia de ADN")
        st.pyplot(fig)

    elif ilustracion == 'Frecuencia de Codones':
        graficar_codones_interactivo(secuencia_adn)

    elif ilustracion == 'Secuencias de Replicación':
        secuencias_replicacion = analizar_replicacion(secuencia_adn)
        if secuencias_replicacion:
            st.write("Se han encontrado las siguientes secuencias de replicación (orígenes de replicación):")
            st.write(secuencias_replicacion)
        else:
            st.write("No se encontraron secuencias de replicación en la secuencia.")

    elif ilustracion == 'Secuencias de Regulación':
        secuencias_reguladoras = detectar_secuencias_reguladoras(secuencia_adn)
        if secuencias_reguladoras:
            st.write("Se han encontrado las siguientes secuencias reguladoras:")
            st.write(secuencias_reguladoras)
        else:
            st.write("No se encontraron secuencias reguladoras en la secuencia.")

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
