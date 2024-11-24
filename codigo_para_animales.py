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

# Función para graficar un análisis de la longitud de secuencias
def graficar_longitud_secuencias():
    longitudes = {animal: len(secuencia) for animal, secuencia in secuencias_adn.items()}

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(longitudes.keys(), longitudes.values(), color='skyblue')
    ax.set_xlabel('Animal')
    ax.set_ylabel('Longitud de Secuencia (nucleótidos)')
    ax.set_title('Longitud de Secuencias de ADN de Diferentes Animales')
    st.pyplot(fig)

# Función para simular mutaciones aleatorias en las secuencias de ADN
def graficar_mutaciones():
    # Simulación de mutaciones: seleccionamos una base al azar y la mutamos
    def mutar_secuencia(secuencia):
        pos = np.random.randint(0, len(secuencia))  # Posición aleatoria en la secuencia
        base_original = secuencia[pos]
        bases = ['A', 'T', 'C', 'G']
        bases.remove(base_original)  # Excluir la base original
        nueva_base = np.random.choice(bases)  # Elegir una nueva base aleatoria
        mutada = secuencia[:pos] + nueva_base + secuencia[pos+1:]
        return mutada

    secuencia_mutada = mutar_secuencia(secuencias_adn['Perro'])  # Ejemplo: mutar la secuencia del Perro

    # Mostrar antes y después de la mutación
    st.write(f"Secuencia Original del Perro: {secuencias_adn['Perro']}")
    st.write(f"Secuencia Mutada: {secuencia_mutada}")
    st.write("Se ha mutado un solo nucleótido aleatorio en la secuencia.")

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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Distribución de Secuencias', 'Longitud de Secuencias', 'Mutaciones Aleatorias']
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

    elif ilustracion == 'Longitud de Secuencias':
        graficar_longitud_secuencias()

    elif ilustracion == 'Mutaciones Aleatorias':
        graficar_mutaciones()

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
