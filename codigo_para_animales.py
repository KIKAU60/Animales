import streamlit as st
from funciones_adn import (
    generar_helice_adn_interactiva,
    calcular_frecuencia_aminos,
    graficar_aminos,
    generar_mapa_genomico,
    generar_grafico_dispersión
)
from secuencias_adn import secuencias_adn

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
    opcion = st.selectbox("Elige una opción para ver:", ["Ver ADN en 3D", "Frecuencia de Aminoácidos", "Mapa Genómico", "Gráfico de Dispersión"])

    if opcion == "Ver ADN en 3D":
        generar_helice_adn_interactiva(secuencia_adn)
    elif opcion == "Frecuencia de Aminoácidos":
        frecuencia_aminos = calcular_frecuencia_aminos(secuencia_adn)
        graficar_aminos(frecuencia_aminos)
    elif opcion == "Mapa Genómico":
        generar_mapa_genomico()
    elif opcion == "Gráfico de Dispersión":
        generar_grafico_dispersión()

# Ejecutar la aplicación
if __name__ == "__main__":
    main()
