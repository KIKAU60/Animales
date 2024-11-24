import sys
import streamlit as st

# Asegurarse de que la ruta al directorio sea la correcta
# Si la ruta está correcta, esto asegura que `secuencias_adn` sea accesible.
sys.path.append('/ruta/al/directorio')

# Importar las funciones desde archivos separados
from funciones_adn import (
    generar_helice_adn_interactiva,
    calcular_frecuencia_aminos,
    graficar_aminos,
    generar_mapa_genomico,
    generar_grafico_dispersión
)

# Importar el diccionario de secuencias de ADN desde el archivo de secuencias
try:
    from secuencias_adn import secuencias_adn
except ImportError as e:
    st.error(f"Error al importar las secuencias de ADN: {e}")
    secuencias_adn = {}  # En caso de que no se importe, asegurar que el diccionario esté vacío

# Función principal de la aplicación Streamlit
def main():
    # Título de la aplicación
    st.title("Análisis de ADN de Animales")

    if not secuencias_adn:
        st.error("No se pudieron cargar las secuencias de ADN. Verifique el archivo 'secuencias_adn.py'.")
        return  # Detener la ejecución si no se pudo cargar el ADN

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
