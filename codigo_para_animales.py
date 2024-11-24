import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
from io import StringIO
import seaborn as sns
import scipy.cluster.hierarchy as sch
from sklearn.metrics import pairwise_distances

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

# Función para crear un heatmap de las secuencias de ADN
def visualizar_secuencias_heatmap():
    # Crear una matriz de distancias de Hamming entre las secuencias
    secuencias = list(secuencias_adn.values())
    distancias = np.zeros((len(secuencias), len(secuencias)))
    
    for i in range(len(secuencias)):
        for j in range(len(secuencias)):
            distancias[i, j] = sum([1 for a, b in zip(secuencias[i], secuencias[j]) if a != b])
    
    # Heatmap con Seaborn
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(distancias, xticklabels=list(secuencias_adn.keys()), yticklabels=list(secuencias_adn.keys()), cmap="YlGnBu", annot=True, fmt="d", ax=ax)
    ax.set_title("Mapa de Distancia entre Secuencias de ADN")
    st.pyplot(fig)

# Función para graficar un mapa de similaridad de secuencias
def mapa_similaridad_secuencias():
    # Convertir las secuencias a formato binario para comparar
    def secuencia_a_binaria(secuencia):
        return [1 if base == 'A' else 0 for base in secuencia]  # Convertir A -> 1, T -> 0
    
    secuencias_binarias = [secuencia_a_binaria(secuencia) for secuencia in secuencias_adn.values()]
    
    # Calcular distancias entre secuencias (utilizando Hamming o alguna métrica binaria)
    distancias = pairwise_distances(secuencias_binarias, metric='hamming')
    
    # Realizar clustering jerárquico
    linkage = sch.linkage(distancias, method='average')
    dendrogram = sch.dendrogram(linkage, labels=list(secuencias_adn.keys()), orientation='top')
    
    st.pyplot()

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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Heatmap de Secuencias', 'Mapa de Similaridad']
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

    elif ilustracion == 'Heatmap de Secuencias':
        visualizar_secuencias_heatmap()

    elif ilustracion == 'Mapa de Similaridad':
        mapa_similaridad_secuencias()

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
