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
    "Elefante": "ATGAGCTGAGAGTCCAGGCGTCGAGGGAGGCGTAGAGGAAGGCGAGT",
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

# Función para analizar enlaces hidrófobos
def analizar_enlaces_hidrofobos(secuencia):
    # Definir regiones hidrófobas como aquellas que contienen "A", "V", "I", "L", "M", "F", "W", "Y"
    regiones_hidrofobas = re.findall(r'[AVILMFYW]{3,}', secuencia)
    return regiones_hidrofobas

# Función para analizar la longitud de las secuencias
def analizar_longitud_secuencias():
    longitudes = [len(secuencia) for secuencia in secuencias_adn.values()]
    promedio_longitud = np.mean(longitudes)
    desviacion_longitud = np.std(longitudes)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.hist(longitudes, bins=10, color='skyblue', edgecolor='black')
    ax.axvline(promedio_longitud, color='red', linestyle='dashed', linewidth=1)
    ax.axvline(promedio_longitud + desviacion_longitud, color='green', linestyle='dashed', linewidth=1)
    ax.axvline(promedio_longitud - desviacion_longitud, color='green', linestyle='dashed', linewidth=1)
    ax.set_title("Distribución de Longitudes de Secuencias de ADN")
    ax.set_xlabel("Longitud de la Secuencia")
    ax.set_ylabel("Frecuencia")
    st.pyplot(fig)

# Función para simular mutaciones genéticas
def simular_mutaciones(secuencia, num_mutaciones=5):
    secuencia_mutada = list(secuencia)
    posiciones = random.sample(range(len(secuencia)), num_mutaciones)
    nucleotidos = ['A', 'T', 'C', 'G']
    
    for pos in posiciones:
        nucleotido_original = secuencia[pos]
        nuevo_nucleotido = random.choice([n for n in nucleotidos if n != nucleotido_original])
        secuencia_mutada[pos] = nuevo_nucleotido
    
    return ''.join(secuencia_mutada)

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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Análisis de Enlaces Hidrófobos', 'Longitud de Secuencias', 'Simulación de Mutaciones']
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

    elif ilustracion == 'Análisis de Enlaces Hidrófobos':
        enlaces_hidrofobos = analizar_enlaces_hidrofobos(secuencia_adn)
        if enlaces_hidrofobos:
            st.write("Se han encontrado las siguientes regiones hidrófobas:")
            st.write(enlaces_hidrofobos)
        else:
            st.write("No se encontraron regiones hidrófobas en la secuencia.")

    elif ilustracion == 'Longitud de Secuencias':
        analizar_longitud_secuencias()

    elif ilustracion == 'Simulación de Mutaciones':
        secuencia_mutada = simular_mutaciones(secuencia_adn)
        st.write("Secuencia Original:")
        st.text(secuencia_adn)
        st.write("Secuencia Tras Mutaciones Aleatorias:")
        st.text(secuencia_mutada)

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
