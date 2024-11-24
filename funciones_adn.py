import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from Bio.Seq import Seq
import pandas as pd
import plotly.graph_objs as go
import plotly.express as px

# Función para graficar la representación 3D de la doble hélice de ADN interactiva
def generar_helice_adn_interactiva(secuencia_adn):
    colores = {'A': 'blue', 'T': 'red', 'C': 'green', 'G': 'yellow'}
    t = np.linspace(0, 4 * np.pi, len(secuencia_adn))  # Usamos 4 pi para 2 vueltas completas
    x = np.sin(t)  # Coordenada X
    y = np.cos(t)  # Coordenada Y
    z = np.linspace(0, 1, len(secuencia_adn))  # Coordenada Z

    # Crear la figura interactiva con plotly
    trace = go.Scatter3d(
        x=x, y=y, z=z,
        mode='markers',
        marker=dict(
            size=5,
            color=[colores[base] for base in secuencia_adn],
            colorscale='Viridis',
            opacity=0.8
        ),
        text=secuencia_adn,
    )

    layout = go.Layout(
        title="Representación 3D de la Doble Hélice de ADN",
        scene=dict(
            xaxis=dict(title="X"),
            yaxis=dict(title="Y"),
            zaxis=dict(title="Z")
        ),
        margin=dict(l=0, r=0, b=0, t=40)
    )

    fig = go.Figure(data=[trace], layout=layout)
    st.plotly_chart(fig, use_container_width=True)  # Mostrar el gráfico interactivo

# Función para calcular la frecuencia de los aminoácidos
def calcular_frecuencia_aminos(secuencia_adn):
    # Convertir ADN a ARN y luego a una secuencia de aminoácidos (simplificado)
    secuencia_arn = Seq(secuencia_adn).transcribe()
    secuencia_aminos = str(secuencia_arn.translate())

    # Contar la frecuencia de los aminoácidos
    frecuencia_aminos = Counter(secuencia_aminos)
    
    return frecuencia_aminos

# Función para graficar los aminoácidos en un gráfico de barras
def graficar_aminos(frecuencia_aminos):
    amino_acidos = list(frecuencia_aminos.keys())
    frecuencias = list(frecuencia_aminos.values())
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.bar(amino_acidos, frecuencias, color='purple')
    ax.set_title('Frecuencia de Aminoácidos')
    ax.set_xlabel('Aminoácido')
    ax.set_ylabel('Frecuencia')
    st.pyplot(fig)  # Mostrar la gráfica con Streamlit

# Función para generar un mapa interactivo de características genómicas
def generar_mapa_genomico():
    # Datos de ejemplo sobre características genómicas (esto debe ser adecuado a tu contexto)
    datos = {
        "Animal": ["Perro", "Gato", "Elefante", "Delfín", "Caballo", "León", "Oso", "Tigre", "Zebra", "Koala"],
        "Latitude": [34.0522, 40.7128, 10.8231, 28.5383, 41.8781, -33.8688, 45.4215, 25.7617, -33.8688, -27.4698],
        "Longitude": [-118.2437, -74.0060, 106.6297, -81.3792, -87.6298, 151.2093, -75.6992, -80.1918, 151.2093, 153.0210]
    }

    df = pd.DataFrame(datos)

    # Usar Plotly Express para crear un mapa interactivo
    fig = px.scatter_geo(df, lat='Latitude', lon='Longitude', text='Animal', 
                         title="Mapa de Características Genómicas de Animales", 
                         template="plotly", scope="world")
    st.plotly_chart(fig)

# Función para crear un gráfico de dispersión
def generar_grafico_dispersión():
    # Datos de ejemplo para el gráfico de dispersión (esto debe ser adecuado a tu contexto)
    x = np.random.rand(100)
    y = np.random.rand(100)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.scatter(x, y, color='blue')
    ax.set_title('Gráfico de Dispersión Genómica')
    ax.set_xlabel('Componente 1')
    ax.set_ylabel('Componente 2')
    st.pyplot(fig)
