import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq

# Diccionario con secuencias de ADN de 10 animales (esto es solo un ejemplo)
# Las secuencias deben ser cadenas de bases A, T, C, G
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
    "Koala": "ATGCGTACGTTAGCCTAGCTAGGCTAGGCTC"
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

# Función para calcular la cantidad de proteínas codificadas
def calcular_proteinas(secuencia_adn):
    secuencia = Seq(secuencia_adn)
    proteina = secuencia.translate()  # Traducción de la secuencia de ADN a proteína
    return proteina

# Función principal de la aplicación Streamlit
def main():
    # Título de la aplicación
    st.title("Visualización de Nucleótidos y Proteínas en el ADN")

    # Crear un selector para elegir entre los 10 animales
    animal = st.selectbox("Selecciona un animal:", list(secuencias_adn.keys()))

    # Obtener la secuencia de ADN del animal seleccionado
    secuencia_adn = secuencias_adn[animal]

    # Mostrar la secuencia de ADN seleccionada
    st.write(f"Secuencia de ADN del {animal}: {secuencia_adn}")

    # Graficar la cantidad de nucleótidos (A, T, C, G)
    graficar_nucleotidos(secuencia_adn)

    # Calcular y mostrar la cantidad de proteínas
    proteina = calcular_proteinas(secuencia_adn)
    st.write(f"La proteína codificada por el ADN del {animal} es: {proteina}")

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
