## TFM_bioinformatica_laurafortaznar
TFM bioinformática por la UNIR

Este repositorio contiene el Trabajo Fin de Máster en Bioinformática. Aquí se incluyen:

- El documento en LaTeX (`main.tex`)
- Los scripts de análisis bioinformático
- Resultados, gráficos y figuras
- Bibliografía en formato BibTeX

## Estructura del proyecto

# Título: Detección de Biomarcadores en Pacientes de Alzheimer Relacionados con Retrovirus Endógenos

Este repositorio contiene el **Trabajo Fin de Máster (TFM)** que tiene como objetivo la **detección de biomarcadores** en pacientes de **Alzheimer** relacionados con los **retrovirus endógenos**. El análisis se realizará utilizando datos de **snRNA-Seq** (RNA secuenciado a nivel de núcleo único) extraídos de astrocitos del córtex de cerebros de pacientes con Alzheimer y sujetos sanos.

---

## Objetivos del TFM

- **Aprender a trabajar con bases de datos en R**.
- **Relacionar diferentes vías moleculares** gracias a la utilización de algoritmos de análisis de datos.
- **Estratificar entre pacientes de Alzheimer genético y esporádico**, y entre sujetos sanos.
- **Utilizar los datos de este artículo**: 
  - *Decoding the Role of Astrocytes in the Entorhinal Cortex in Alzheimer's Disease Using High-Dimensional Single-Nucleus RNA Sequencing Data and Next-Generation Knowledge Discovery Methodologies: Focus on Drugs and Natural Product Remedies for Dementia*  
  - **Autores**: Peter Natesan Pushparaj, Gauthaman Kalamegam, Khalid Hussain Wali Sait, Mahmood Rasool  
  - **PMID**: [35295737](https://pubmed.ncbi.nlm.nih.gov/35295737/)  
  - **PMCID**: [PMC8918735](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8918735/)  
  - **DOI**: [10.3389/fphar.2021.720170](https://doi.org/10.3389/fphar.2021.720170)

---

## 📁 Estructura del repositorio

📂 TFM-Bioinformatica/ ├── main.tex # Documento principal en LaTeX ├── references.bib # Archivo de bibliografía en BibTeX ├── .gitignore # Archivos ignorados por Git ├── README.md # Este archivo │ ├── sections/ # Capítulos separados │ ├── intro.tex │ ├── methods.tex │ ├── results.tex │ └── discussion.tex │ ├── images/ # Figuras y gráficos generados durante el análisis │ ├── figura1.png │ └── ... │ ├── scripts/ # Scripts de análisis (R, Python, etc.) │ ├── analisis.R │ └── estratificacion.R │ └── data/ # Datos utilizados en el análisis ├── raw/ # Datos crudos del artículo └── processed/ # Datos procesados para análisis


## 🛠️ Tecnologías utilizadas

- **R** para el análisis de datos y bioinformática
- **LaTeX** para la redacción científica
- **Git** + **GitHub** para control de versiones
- **snRNA-Seq** para datos de expresión génica a nivel de núcleo único
- Bibliografía gestionada con **BibTeX**

---

## 📌 Notas

- Este proyecto usa datos de los artículos mencionados y de la base de datos **GEO (GSE138852, GSE147528)**.
- La compilación del documento LaTeX debe realizarse con `pdflatex`, `xelatex`, o `latexmk` para obtener el archivo PDF final.
- Los scripts están en R y están organizados para realizar análisis específicos relacionados con la estratificación de los pacientes y biomarcadores.

---

## ✍️ Autor

Laura Fort Aznar – [laurasfort@gmail.com]  
Máster en Bioinformática  
[UNIR, 2025]

---

## 📝 Licencia

Este proyecto está bajo la licencia pública general GNU v3 (GPLv3).
