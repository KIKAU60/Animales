import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
import seaborn as sns
from io import StringIO
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram, linkage

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

# Función para analizar la distribución de las secuencias
def graficar_distribucion():
    # Calcular proporciones de A, T, C, G para cada secuencia
    proporciones = {}
    for animal, secuencia in secuencias_adn.items():
        proporciones[animal] = calcular_proporcion_nucleotidos(secuencia)
    
    # Convertir en formato adecuado para gráfico
    labels = ['A', 'T', 'C', 'G']
    distribucion = np.array([list(proporciones[animal]) for animal in secuencias_adn])

    # Graficar la distribución de bases con un gráfico de dispersión
    fig, ax = plt.subplots(figsize=(8, 8))
    scatter = ax.scatter(distribucion[:, 0], distribucion[:, 1], c=distribucion[:, 2], s=100, cmap='viridis')
    ax.set_xlabel('Proporción A')
    ax.set_ylabel('Proporción T')
    ax.set_title('Distribución de Secuencias de ADN por Proporción de Bases')
    fig.colorbar(scatter, ax=ax, label='Proporción C')
    ax.legend([animal for animal in secuencias_adn], title='Animales')
    st.pyplot(fig)

# Función para graficar un árbol filogenético
def graficar_arbol_filogenetico():
    # Crear una matriz de distancias de Hamming entre las secuencias
    secuencias = list(secuencias_adn.values())
    distancias = np.zeros((len(secuencias), len(secuencias)))
    
    for i in range(len(secuencias)):
        for j in range(len(secuencias)):
            distancias[i, j] = sum([1 for a, b in zip(secuencias[i], secuencias[j]) if a != b])

    # Usar linkage de scipy para calcular la jerarquía
    linked = linkage(distancias, method='average')
    
    # Graficar el dendrograma
    fig, ax = plt.subplots(figsize=(10, 8))
    dendrogram(linked, labels=list(secuencias_adn.keys()), orientation='top', distance_sort='descending', show_leaf_counts=True)
    ax.set_title("Árbol Filogenético de Secuencias de ADN")
    st.pyplot(fig)

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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Distribución de Secuencias', 'Árbol Filogenético']
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

    elif ilustracion == 'Distribución de Secuencias':
        graficar_distribucion()

    elif ilustracion == 'Árbol Filogenético':
        graficar_arbol_filogenetico()

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
