import sys
import streamlit as st

# Añadir la ruta al directorio donde se encuentra 'secuencias_adn.py'
sys.path.append('')

# Ahora puedes importar el archivo 'secuencias_adn' correctamente
from secuencias_de_adn import secuencias_de_adn

# Importar las funciones necesarias
from funciones_adn import (
    generar_helice_adn_interactiva,
    calcular_frecuencia_aminos,
    graficar_aminos,
    generar_mapa_genomico,
    generar_grafico_dispersión
)

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
