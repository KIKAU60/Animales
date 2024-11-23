import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from Bio.Seq import Seq
from Bio import Entrez
from collections import Counter
import py3Dmol
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import seaborn as sns

# Diccionario con secuencias de ADN de 15 animales
secuencias_adn = {
    "Perro": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTA",
    "Gato": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTT",
    "Elefante": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCGT",
    "Delfín": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCAA",
    "Caballo": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTG",
    "León": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTC",
    "Oso": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTG",
    "Tigre": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTT",
    "Zebra": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTG",
    "Koala": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTC",
    "Llama": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTA",
    "Panda": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTG",
    "Canguro": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTA",
    "Gorila": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTT",
    "León Marino": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTC"
}

# Función para graficar la cantidad de nucleótidos
def graficar_nucleotidos(secuencia_adn):
    secuencia = Seq(secuencia_adn)
    count_a = secuencia.count('A')
    count_t = secuencia.count('T')
    count_c = secuencia.count('C')
    count_g = secuencia.count('G')

    nucleotidos = ['Adenina (A)', 'Timina (T)', 'Citosina (C)', 'Guanina (G)']
    cantidades = [count_a, count_t, count_c, count_g]
    colores = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.bar(nucleotidos, cantidades, color=colores)
    ax.set_title("Cantidad de Nucleótidos en la Secuencia de ADN")
    ax.set_xlabel('Nucleótidos')
    ax.set_ylabel('Cantidad')
    st.pyplot(fig)  # Mostrar la gráfica con Streamlit

# Función para obtener la transcripción de ARN
def obtener_transcripcion_arn(secuencia_adn):
    secuencia_arn = Seq(secuencia_adn).transcribe()
    return str(secuencia_arn)

# Función para obtener la traducción de ARN a proteínas
def obtener_traduccion_arn(secuencia_arn):
    proteina = Seq(secuencia_arn).translate()
    return str(proteina)

# Función para obtener y graficar la frecuencia de codones
def obtener_codones(secuencia_adn):
    return [secuencia_adn[i:i+3] for i in range(0, len(secuencia_adn), 3)]

def graficar_codones(codones):
    frecuencia_codones = Counter(codones)
    codones, frecuencias = zip(*sorted(frecuencia_codones.items()))
    
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(codones, frecuencias, color='purple')
    ax.set_title('Frecuencia de Codones')
    ax.set_xlabel('Codón')
    ax.set_ylabel('Frecuencia')
    plt.xticks(rotation=90)
    st.pyplot(fig)  # Mostrar la gráfica con Streamlit

# Función para generar la representación 3D de la doble hélice de ADN
def generar_helice_adn(secuencia_adn):
    colores = {'A': 'blue', 'T': 'red', 'C': 'green', 'G': 'yellow'}
    t = np.linspace(0, 4 * np.pi, len(secuencia_adn))  # Usamos 4 pi para 2 vueltas completas
    x = np.sin(t)  # Coordenada X
    y = np.cos(t)  # Coordenada Y
    z = np.linspace(0, 1, len(secuencia_adn))  # Coordenada Z

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Dibuja los puntos de la hélice para cada base
    for i, base in enumerate(secuencia_adn):
        ax.scatter(x[i], y[i], z[i], color=colores[base], s=100, label=base if i == 0 else "")

    ax.set_title("Representación 3D de la Doble Hélice de ADN")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(False)
    ax.view_init(30, 60)  # Ajuste del ángulo de vista

    st.pyplot(fig)  # Mostrar la gráfica en Streamlit

# Función para construir el árbol filogenético basado en distancias de ADN
def graficar_arbol_filogenetico(secuencias):
    distancias = np.array([[sum(1 for a, b in zip(seq1, seq2) if a != b) for seq2 in secuencias] for seq1 in secuencias])
    linked = linkage(distancias, 'single')

    fig, ax = plt.subplots(figsize=(10, 8))
    dendrogram(linked, labels=list(secuencias.keys()))
    ax.set_title("Árbol Filogenético de Animales Basado en ADN")
    st.pyplot(fig)  # Mostrar la gráfica con Streamlit

# Función principal de la aplicación Streamlit
def main():
    st.title("Análisis de ADN de 15 Animales")

    # Crear un selector para elegir entre los 15 animales
    animal = st.selectbox("Selecciona un animal:", list(secuencias_adn.keys()))

    # Obtener la secuencia de ADN del animal seleccionado
    secuencia_adn = secuencias_adn[animal]

    # Mostrar la secuencia de ADN seleccionada
    st.write(f"Secuencia de ADN del {animal}: {secuencia_adn}")

    # Visualizar la doble hélice
    generar_helice_adn(secuencia_adn)

    # Graficar la cantidad de nucleótidos (A, T, C, G)
    graficar_nucleotidos(secuencia_adn)

    # Obtener y graficar la frecuencia de codones
    codones = obtener_codones(secuencia_adn)
    graficar_codones(codones)

    # Calcular y mostrar la transcripción del ARN
    secuencia_arn = obtener_transcripcion_arn(secuencia_adn)
    st.write(f"Secuencia de ARN transcrita: {secuencia_arn}")

    # Calcular y mostrar la traducción a proteína
    secuencia_proteina = obtener_traduccion_arn(secuencia_arn)
    st.write(f"Secuencia de proteína traducida: {secuencia_proteina}")

    # Graficar el árbol filogenético
    graficar_arbol_filogenetico(secuencias_adn)

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
