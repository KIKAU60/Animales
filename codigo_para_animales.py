import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
from io import StringIO
import seaborn as sns
import scipy.cluster.hierarchy as sch

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

# Función para ilustrar la relación filogenética (Árbol)
def ilustrar_evolucion_filogenetica():
    # Crear datos de ejemplo (distancias de Hamming entre secuencias de ADN)
    animals = list(secuencias_adn.keys())
    secuencias = list(secuencias_adn.values())
    
    # Crear una matriz de distancias de Hamming (simuladas en este caso)
    distancias = np.zeros((len(secuencias), len(secuencias)))
    
    for i in range(len(secuencias)):
        for j in range(len(secuencias)):
            distancias[i, j] = sum([1 for a, b in zip(secuencias[i], secuencias[j]) if a != b])
    
    # Crear un árbol jerárquico usando scipy
    linkage = sch.linkage(distancias, method='average')
    dendrogram = sch.dendrogram(linkage, labels=animals, orientation='right')
    
    st.pyplot()

# Función para identificar y graficar motivos conservados
def buscar_motivos_conservados(secuencia, longitud_motivo=5):
    # Buscar motivos conservados de longitud especificada
    motivos = [secuencia[i:i+longitud_motivo] for i in range(len(secuencia) - longitud_motivo + 1)]
    counter = Counter(motivos)
    
    # Mostrar los motivos más frecuentes
    motivos_frecuentes = counter.most_common(10)
    
    # Crear un gráfico de barras con los motivos más frecuentes
    motivos_unicos = [motivo for motivo, _ in motivos_frecuentes]
    frecuencias = [freq for _, freq in motivos_frecuentes]
    
    fig = plt.figure(figsize=(10, 6))
    plt.bar(motivos_unicos, frecuencias, color='lightcoral')
    plt.xticks(rotation=45, ha='right')
    plt.xlabel('Motivo Conservado')
    plt.ylabel('Frecuencia')
    plt.title(f"Motivos Conservados en la Secuencia de ADN")
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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Evolución Filogenética', 'Motivos Conservados']
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

    elif ilustracion == 'Evolución Filogenética':
        ilustrar_evolucion_filogenetica()

    elif ilustracion == 'Motivos Conservados':
        buscar_motivos_conservados(secuencia_adn)

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
