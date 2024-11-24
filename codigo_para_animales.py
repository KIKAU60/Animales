import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
import seaborn as sns
from io import StringIO
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

# Función para analizar secuencias conservadas entre diferentes animales
def buscar_secuencias_conservadas():
    secuencias = list(secuencias_adn.values())
    secuencias_conservadas = []

    # Buscar subsecuencias comunes entre todas las secuencias
    min_longitud = min(len(seq) for seq in secuencias)
    for i in range(min_longitud - 2):
        subsecuencia = secuencias[0][i:i+3]  # Usamos subsecuencias de longitud 3
        if all(subsecuencia in seq for seq in secuencias):
            secuencias_conservadas.append(subsecuencia)
    
    return secuencias_conservadas

# Función para identificar marcos de lectura abiertos (ORFs)
def identificar_orfs(secuencia):
    orfs = []
    for i in range(len(secuencia) - 3):
        if secuencia[i:i+3] == "ATG":  # Iniciar un ORF (codón de inicio)
            for j in range(i+3, len(secuencia)-3, 3):
                if secuencia[j:j+3] == "TAA" or secuencia[j:j+3] == "TAG" or secuencia[j:j+3] == "TGA":  # Codones de parada
                    orfs.append(secuencia[i:j+3])
                    break
    return orfs

# Función para simular la duplicación de ADN a lo largo del tiempo
def duplicacion_adn(secuencia, ciclos=5):
    secuencia_duplicada = secuencia
    for _ in range(ciclos):
        secuencia_duplicada += secuencia
    return secuencia_duplicada

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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Secuencias Conservadas', 'Identificación de ORFs', 'Duplicación de Secuencia']
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

    elif ilustracion == 'Secuencias Conservadas':
        secuencias_conservadas = buscar_secuencias_conservadas()
        if secuencias_conservadas:
            st.write("Las siguientes secuencias están conservadas en todas las especies:")
            st.write(secuencias_conservadas)
        else:
            st.write("No se encontraron secuencias conservadas en todas las especies.")

    elif ilustracion == 'Identificación de ORFs':
        orfs = identificar_orfs(secuencia_adn)
        if orfs:
            st.write("Se han identificado los siguientes marcos de lectura abiertos (ORFs):")
            st.write(orfs)
        else:
            st.write("No se identificaron ORFs en la secuencia.")

    elif ilustracion == 'Duplicación de Secuencia':
        secuencia_duplicada = duplicacion_adn(secuencia_adn)
        st.write("Secuencia Original:")
        st.text(secuencia_adn)
        st.write("Secuencia Duplicada (tras varias copias):")
        st.text(secuencia_duplicada)

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
