import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from collections import Counter
import seaborn as sns
from io import StringIO
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import dendrogram, linkage
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

# Función para buscar palíndromos en la secuencia de ADN
def buscar_palindromos(secuencia):
    # Definir el patrón de palíndromos (secuencias que se leen igual de izquierda a derecha que de derecha a izquierda)
    palindromos = []
    for i in range(len(secuencia)):
        for j in range(i+4, len(secuencia)+1):  # Palíndromos de longitud mínima 4
            sub_secuencia = secuencia[i:j]
            if sub_secuencia == sub_secuencia[::-1]:  # Verifica si es un palíndromo
                palindromos.append(sub_secuencia)
    return palindromos

# Función para replicar una secuencia de ADN (simulación)
def replicar_secuencia(secuencia):
    # Simula la replicación de ADN: simplemente se duplica la secuencia
    return secuencia + secuencia

# Función para detectar motivos de unión de proteínas
def detectar_motivos_union(secuencia):
    # Buscar patrones comunes en la secuencia que podrían ser sitios de unión de proteínas (ejemplo simple)
    # Por ejemplo, buscamos secuencias de 6 bases (que podrían representar motivos de unión)
    patrones = re.findall(r'(?=(ATG[A-Z]{3}G))', secuencia)  # Simple búsqueda de patrones
    return patrones

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
        ['Proporciones de Nucleótidos', 'Frecuencia de Codones', 'Palíndromos en ADN', 'Replicación de Secuencia', 'Motivos de Unión de Proteínas']
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

    elif ilustracion == 'Palíndromos en ADN':
        palindromos = buscar_palindromos(secuencia_adn)
        if palindromos:
            st.write("Se han encontrado los siguientes palíndromos en la secuencia:")
            st.write(palindromos)
        else:
            st.write("No se han encontrado palíndromos en la secuencia.")

    elif ilustracion == 'Replicación de Secuencia':
        secuencia_replicada = replicar_secuencia(secuencia_adn)
        st.write("Secuencia Original:")
        st.text(secuencia_adn)
        st.write("Secuencia Replicada (duplicada):")
        st.text(secuencia_replicada)

    elif ilustracion == 'Motivos de Unión de Proteínas':
        motivos = detectar_motivos_union(secuencia_adn)
        if motivos:
            st.write("Se han encontrado los siguientes motivos de unión de proteínas:")
            st.write(motivos)
        else:
            st.write("No se han encontrado motivos de unión de proteínas en la secuencia.")

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
