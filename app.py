import streamlit as st
import random
import string
import matplotlib.pyplot as plt

# Función para generar una secuencia de ADN aleatoria de un tamaño dado
def generate_random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

# Función para calcular el contenido GC de la secuencia
def calculate_gc_content(sequence):
    gc_count = sum(1 for base in sequence if base in 'GCgc')
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

# Función para analizar y mostrar la secuencia de ADN
def visualize_dna(sequence):
    if sequence:
        sequence = sequence.strip().upper()  # Limpiar y convertir a mayúsculas
        length = len(sequence)
        gc_content = calculate_gc_content(sequence)

        # Mostrar la secuencia de ADN
        st.subheader("Secuencia de ADN:")
        st.text(sequence)

        # Mostrar información adicional
        st.subheader("Información adicional:")
        st.write(f"Longitud de la secuencia: {length} bases")
        st.write(f"Contenido GC: {gc_content:.2f}%")
        st.write("Posibles características genéticas:")
        st.write("- Aumento en la producción de proteínas (si es un gen codificador).")
        st.write("- Variante genética de adaptación a diferentes ambientes.")

        # Mostrar gráfico de contenido GC
        st.subheader("Gráfico del contenido GC")
        fig, ax = plt.subplots()
        ax.bar(["A", "T", "G", "C"], [sequence.count('A'), sequence.count('T'), sequence.count('G'), sequence.count('C')])
        ax.set_ylabel("Número de bases")
        ax.set_title("Distribución de bases en la secuencia de ADN")
        st.pyplot(fig)

    else:
        st.warning("Por favor, introduce una secuencia de ADN válida.")

# Título de la aplicación
st.title("Análisis de Secuencias de ADN Animal")

# Opción para generar una secuencia de ADN aleatoria
generate_sequence = st.checkbox("Generar secuencia de ADN aleatoria")

# Si el usuario quiere generar una secuencia aleatoria
if generate_sequence:
    sequence_length = st.slider("Selecciona la longitud de la secuencia", 10, 1000, 100)
    random_sequence = generate_random_dna(sequence_length)
    st.text_area("Secuencia generada aleatoriamente:", random_sequence, height=150)

# Opción para que el usuario ingrese su propia secuencia de ADN
sequence_input = st.text_area(
    "Introduce una secuencia de ADN:", 
    help="Introduce la secuencia de ADN en formato estándar (ejemplo: ATGCGTACGTTAGC...)."
)

# Botón para visualizar la secuencia
if st.button("Visualizar Secuencia"):
    if generate_sequence:
        visualize_dna(random_sequence)
    else:
        visualize_dna(sequence_input)
