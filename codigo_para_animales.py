import random
import streamlit as st

# Función para generar una secuencia de ADN de longitud n
def generar_adn(longitud=100):
    secuencia = ''.join(random.choices('ACGT', k=longitud))
    return secuencia

# Definir los 4 indicadores, que son secuencias características
indicadores = ["ATGC", "CGTA", "TGCATG", "GCGT"]

# Generar 10 secuencias de ADN con 4 indicadores en cada una
def generar_secuencias_con_indicadores(num_secuencias=10, longitud=100):
    secuencias = []
    for _ in range(num_secuencias):
        secuencia = generar_adn(longitud)
        
        # Inserta aleatoriamente los 4 indicadores en diferentes posiciones
        for indicador in indicadores:
            posicion = random.randint(0, longitud - len(indicador))
            secuencia = secuencia[:posicion] + indicador + secuencia[posicion + len(indicador):]
        
        secuencias.append(secuencia)
    
    return secuencias

# Función para mostrar el ADN con los indicadores
def mostrar_secuencias(secuencias):
    for idx, secuencia in enumerate(secuencias, 1):
        st.subheader(f"Secuencia {idx}: {secuencia}")
        st.write("Indicadores en la secuencia:")
        
        for indicador in indicadores:
            if indicador in secuencia:
                st.write(f"  - {indicador} encontrado.")
        
        st.markdown("---")  # Línea separadora

# Título y Descripción de la App
st.title("Generador de Secuencias de ADN con Indicadores")
st.write("""
    Esta aplicación genera secuencias aleatorias de ADN con indicadores específicos en ellas. 
    Los indicadores son secuencias cortas de ADN que se insertan de manera aleatoria.
""")

# Parámetros para generar las secuencias
longitud_adn = st.slider("Longitud de las secuencias de ADN", min_value=50, max_value=200, value=100)
num_secuencias = st.slider("Número de secuencias a generar", min_value=1, max_value=20, value=10)

# Botón para generar las secuencias
if st.button("Generar Secuencias"):
    secuencias_con_indicadores = generar_secuencias_con_indicadores(num_secuencias=num_secuencias, longitud=longitud_adn)
    mostrar_secuencias(secuencias_con_indicadores)

