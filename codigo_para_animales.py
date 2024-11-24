import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter

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

# Función para ilustrar la estructura de la doble hélice (gráfico 3D simple)
def ilustrar_doble_helice(secuencia):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Simulación de puntos en 3D (esto es solo una representación visual)
    x = np.random.rand(10)
    y = np.random.rand(10)
    z = np.random.rand(10)
    
    ax.plot(x, y, z, label="Estructura de la Doble Hélice")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title("Estructura de la Doble Hélice de ADN")
    ax.legend()
    st.pyplot(fig)

# Función para mostrar la secuencia de ADN como texto
def representar_secuencia(secuencia):
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.text(0.5, 0.5, secuencia, fontsize=12, ha='center', va='center')
    ax.set_title("Representación de la Secuencia de ADN")
    ax.axis('off')
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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Estructura de Doble Hélice', 'Representación de la Secuencia']
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

    elif ilustracion == 'Estructura de Doble Hélice':
        ilustrar_doble_helice(secuencia_adn)

    elif ilustracion == 'Representación de la Secuencia':
        representar_secuencia(secuencia_adn)

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
