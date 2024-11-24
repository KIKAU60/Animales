import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import Counter
from Bio.Seq import Seq
import pandas as pd
from matplotlib.animation import FuncAnimation

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

# Función para graficar la representación 3D de la doble hélice de ADN con animación
def generar_helice_adn_animada(secuencia_adn):
    colores = {'A': 'blue', 'T': 'red', 'C': 'green', 'G': 'yellow'}
    t = np.linspace(0, 4 * np.pi, len(secuencia_adn))  # Usamos 4 pi para 2 vueltas completas
    x = np.sin(t)  # Coordenada X
    y = np.cos(t)  # Coordenada Y
    z = np.linspace(0, 1, len(secuencia_adn))  # Coordenada Z

    fig = plt.figure(figsize=(10, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Dibuja los puntos de la hélice para cada base
    scatters = []
    for i, base in enumerate(secuencia_adn):
        scatter = ax.scatter(x[i], y[i], z[i], color=colores[base], s=100)
        scatters.append(scatter)

    # Conectar las bases complementarias entre las dos cadenas
    for i in range(0, len(secuencia_adn) - 1, 2):
        ax.plot([x[i], x[i+1]], [y[i], y[i+1]], [z[i], z[i+1]], color='black', lw=1)

    ax.set_title("Representación 3D de la Doble Hélice de ADN")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid(False)
    ax.view_init(30, 60)  # Ajuste inicial del ángulo de vista

    # Función de actualización para la animación (giro)
    def actualizar(i):
        ax.view_init(elev=30, azim=i)  # Cambia el ángulo de vista para crear el giro

    # Crear la animación
    ani = FuncAnimation(fig, actualizar, frames=np.arange(0, 360, 2), interval=50)

    # Mostrar la animación en Streamlit
    st.pyplot(fig)  # Mostrar la figura estática por defecto (para asegurar compatibilidad)

# Función principal de la aplicación Streamlit
def main():
    # Título de la aplicación
    st.title("Análisis de ADN de Animales")

    # Crear un selector para elegir entre los 15 animales
    animal = st.selectbox("Selecciona un animal:", list(secuencias_adn.keys()))

    # Obtener la secuencia de ADN del animal seleccionado
    secuencia_adn = secuencias_adn[animal]

    # Mostrar la secuencia de ADN seleccionada
    st.write(f"Secuencia de ADN del {animal}: {secuencia_adn}")

    # Mostrar la representación 3D de la doble hélice con animación de giro
    generar_helice_adn_animada(secuencia_adn)

    # Calcular y mostrar la proporción de nucleótidos
    calcular_proporcion_nucleotidos(secuencia_adn)

    # Obtener y mostrar los codones
    codones = obtener_codones(secuencia_adn)
    graficar_codones(codones)

    # Mostrar la tabla de codones
    mostrar_tabla_codones(codones)

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
