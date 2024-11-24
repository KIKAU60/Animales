import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from Bio.Seq import Seq
import pandas as pd

# Diccionario con secuencias de ADN de 15 animales (solo ejemplo)
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
    "Cebra": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTT",
    "Liebre": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCCA",
    "Mono": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCAG",
    "Cabra": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCGC",
    "Rata": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTA"
}

# Función para graficar la representación 3D de la doble hélice de ADN
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

    # Conectar las bases complementarias entre las dos cadenas
    for i in range(0, len(secuencia_adn) - 1, 2):
        ax.plot([x[i], x[i+1]], [y[i], y[i+1]], [z[i], z[i+1]], color='black', lw=1)

    ax.set_title("Representación 3D de la Doble Hélice de ADN")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(False)
    ax.view_init(30, 60)  # Ajuste del ángulo de vista

    st.pyplot(fig)  # Mostrar la gráfica con Streamlit

# Función para calcular la proporción de nucleótidos (A, T, C, G)
def calcular_proporcion_nucleotidos(secuencia_adn):
    secuencia = Seq(secuencia_adn)
    count_a = secuencia.count('A')
    count_t = secuencia.count('T')
    count_c = secuencia.count('C')
    count_g = secuencia.count('G')
    total_nucleotidos = len(secuencia)

    proporcion_a = count_a / total_nucleotidos
    proporcion_t = count_t / total_nucleotidos
    proporcion_c = count_c / total_nucleotidos
    proporcion_g = count_g / total_nucleotidos

    nucleotidos = ['Adenina (A)', 'Timina (T)', 'Citosina (C)', 'Guanina (G)']
    proporciones = [proporcion_a, proporcion_t, proporcion_c, proporcion_g]
    colores = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.pie(proporciones, labels=nucleotidos, autopct='%1.1f%%', startangle=140, colors=colores)
    ax.set_title("Proporción de Nucleótidos en la Secuencia de ADN")
    st.pyplot(fig)  # Mostrar la gráfica con Streamlit

# Función para obtener los codones de una secuencia de ADN
def obtener_codones(secuencia_adn):
    return [secuencia_adn[i:i+3] for i in range(0, len(secuencia_adn), 3)]

# Función para graficar los codones
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

# Función para crear una tabla con los codones
def mostrar_tabla_codones(codones):
    frecuencia_codones = Counter(codones)
    codones, frecuencias = zip(*sorted(frecuencia_codones.items()))
    df = pd.DataFrame(list(zip(codones, frecuencias)), columns=["Codón", "Frecuencia"])
    st.write(df)

# Función principal de la aplicación
