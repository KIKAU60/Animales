# funciones_adn.py

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# Función para generar una hélice de ADN interactiva (esto sería un ejemplo simple)
def generar_helice_adn_interactiva(secuencia_adn):
    st.write("Generando hélice 3D de ADN...")
    # Aquí usarías una librería de visualización 3D, como py3Dmol o similar
    st.text("Visualización de la hélice de ADN. (Aquí debería ir la visualización interactiva)")

# Función para calcular la frecuencia de aminoácidos en la secuencia
def calcular_frecuencia_aminos(secuencia_adn):
    # Supongamos que la secuencia se traduce a proteínas (simplificación)
    aminoacidos = {'A': 0, 'T': 0, 'C': 0, 'G': 0}  # Aquí se podría hacer una traducción real a aminoácidos
    for base in secuencia_adn:
        if base in aminoacidos:
            aminoacidos[base] += 1
    return aminoacidos

# Función para graficar la frecuencia de aminoácidos
def graficar_aminos(frecuencia_aminos):
    aminoacidos = list(frecuencia_aminos.keys())
    frecuencias = list(frecuencia_aminos.values())
    
    fig, ax = plt.subplots()
    ax.bar(aminoacidos, frecuencias)
    ax.set_xlabel('Aminoácidos')
    ax.set_ylabel('Frecuencia')
    ax.set_title('Frecuencia de Aminoácidos')
    
    st.pyplot(fig)

# Función para generar un mapa genómico (simplificado)
def generar_mapa_genomico():
    st.write("Generando mapa genómico...")
    # Este es solo un ejemplo, el mapa genómico real podría ser interactivo
    st.text("Aquí debería ir un mapa genómico interactivo")

# Función para generar un gráfico de dispersión (simplificado)
def generar_grafico_dispersión():
    st.write("Generando gráfico de dispersión...")
    # Este es solo un ejemplo con datos aleatorios
    x = np.random.rand(50)
    y = np.random.rand(50)
    
    fig, ax = plt.subplots()
    ax.scatter(x, y)
    ax.set_xlabel('Eje X')
    ax.set_ylabel('Eje Y')
    ax.set_title('Gráfico de Dispersión')
    
    st.pyplot(fig)

