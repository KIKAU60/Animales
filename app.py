import streamlit as st
from Bio import Entrez
from Bio import SeqIO

# Configuración de Entrez
Entrez.email = "a223211679@unison.mx"  # Reemplaza con tu correo electrónico

# Función para obtener la secuencia de ADN desde GenBank
def obtener_secuencia_genbank(accession_number):
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        secuencia_adn = str(record.seq)
        return secuencia_adn
    except Exception as e:
        st.error(f"Error al obtener la secuencia de GenBank: {e}")
        return None

# Función principal para la aplicación Streamlit
def main():
    st.title("Obtención de Secuencias de ADN desde GenBank")
    
    # Entrada para el número de acceso de GenBank
    accession_number = st.text_input("Introduce el número de acceso de GenBank (Accession Number):")
    
    if accession_number:
        # Obtener la secuencia desde GenBank
        secuencia_adn = obtener_secuencia_genbank(accession_number)
        
        if secuencia_adn:
            st.write(f"Secuencia de ADN obtenida de GenBank (Accession Number: {accession_number}):")
            st.text(secuencia_adn)
        else:
            st.write("No se pudo obtener la secuencia. Verifica el número de acceso.")

# Ejecutar la aplicación
if __name__ == "__main__":
    main()

