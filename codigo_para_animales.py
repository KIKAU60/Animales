import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from mpl_toolkits.mplot3d import Axes3D

# Secuencias de ADN de animales
secuencias_adn = {
    "Elefante": "AGCTGACGTAGCGTACGTAAGCTG...",
    "Perro": "ATCGAGCTGGTAGCGGATCGAAGT...",
    "Gato": "AAGGCTAGCTAGGTACGTCGAAGTC...",
    "Caballo": "AGGTCGACGTTGAGTCTGAGTGA...",
    "León": "ATGCGATCGTACGAGTGTAGCTAG...",
    "Tigre": "CTGAGTGAGTCGATAGCGATGCAG...",
    "Delfín": "AGTCTGATCGGAGTCTACGAGAGT...",
    "Ballena": "ACGTGAGTACGAGTGTACGTAGTG...",
    "Cebra": "ATGAGTCTAGGATCGAGTACGAGGT...",
    "Rinoceronte": "AGTCGTAGGCTAGCTGACGTAGCG..."
}

# Función para generar la hélice de ADN interactiva (simulada)
def generar_helice_adn_interactiva(secuencia_adn):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(np.random.rand(10), np.random.rand(10), np.random.rand(10), label="ADN en 3D")
    ax.set_title("Estructura de la Doble Hélice de ADN en 3D")
    ax.legend()
    st.pyplot(fig)

# Función para calcular la frecuencia de aminoácidos en la secuencia de ADN
def calcular_frecuencia_aminos(secuencia_adn):
    # Ejemplo de cálculo de frecuencias de aminoácidos
    aminos = ['A', 'T', 'C', 'G']
    frecuencias = {amino: secuencia_adn.count(amino) for amino in aminos}
    return frecuencias

# Función para graficar la frecuencia de aminoácidos
def graficar_aminos(frecuencia_aminos):
    amino_acidos = list(frecuencia_aminos.keys())
    valores = list(frecuencia_aminos.values())
    
    plt.bar(amino_acidos, valores, color='skyblue')
    plt.xlabel('Aminoácidos')
    plt.ylabel('Frecuencia')
    plt.title('Frecuencia de Aminoácidos')
    st.pyplot(plt)

# Función para generar un mapa genómico (gráfico simulado)
def generar_mapa_genomico():
    sns.set(style="whitegrid")
    data = np.random.rand(10, 10)
    sns.heatmap(data, annot=True, cmap="coolwarm")
    plt.title("Mapa Genómico")
    st.pyplot(plt)

# Función para generar un gráfico de dispersión de las secuencias
def generar_grafico_dispersión():
    x = np.random.rand(50)
    y = np.random.rand(50)
    plt.scatter(x, y, color='red')
    plt.xlabel('Coordenada X')
    plt.ylabel('Coordenada Y')
    plt.title('Gráfico de Dispersión Genómica')
    st.pyplot(plt)

# Función principal de la aplicación Streamlit
def main():
    # Título de la aplicación
    st.title("Análisis de ADN de Animales")

    # Crear un selector para elegir entre los 10 animales
    animal = st.selectbox("Selecciona un animal:", list(secuencias_adn.keys()))

    # Obtener la secuencia de ADN del animal seleccionado
    secuencia_adn = secuencias_adn[animal]

    # Mostrar la secuencia de ADN seleccionada
    st.write(f"Secuencia de ADN del {animal}: {secuencia_adn}")

    # Agregar indicadores (opciones) para mostrar diferentes gráficos
    opcion = st.selectbox(
        "Elige una opción para ver:",
        ["Ver ADN en 3D", "Frecuencia de Aminoácidos", "Mapa Genómico", "Gráfico de Dispersión"]
    )

    # Ejecutar la opción seleccionada
    if opcion == "Ver ADN en 3D":
        st.subheader("Visualización de la Doble Hélice de ADN en 3D")
        generar_helice_adn_interactiva(secuencia_adn)
    elif opcion == "Frecuencia de Aminoácidos":
        st.subheader("Frecuencia de los Aminoácidos derivados del ADN")
        frecuencia_aminos = calcular_frecuencia_aminos(secuencia_adn)
        graficar_aminos(frecuencia_aminos)
    elif opcion == "Mapa Genómico":
        st.subheader("Mapa Genómico Interactivo")
        generar_mapa_genomico()
    elif opcion == "Gráfico de Dispersión":
        st.subheader("Gráfico de Dispersión Genómica")
        generar_grafico_dispersión()

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
