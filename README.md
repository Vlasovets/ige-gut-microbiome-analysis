
# The role of IgE patterns and their link to the gut microbiome in allergic sensitization

<p align="center">
	<img src="figures/abstract.png" alt="Visual Abstract" width="500"/>
</p>



This repository contains the analysis code, notebooks, and figures for the study of IgE patterns and their association with the gut microbiome in allergic sensitization.


## Repository Structure

- `notebooks/`: Jupyter notebooks and R Markdown for all analyses
	- `01_IgE_analysis.ipynb`: IgE data analysis and workflow
	- `02_diversity.ipynb`: Microbiome diversity analysis
	- `03_DA_analysis.ipynb`: Differential abundance analysis
	- `04_classification.ipynb`: Classification modeling
	- `05_Network_analysis.ipynb`: Microbial network analysis
	- `06_Network_visualization.Rmd`: Network visualization (R Markdown)
- `figures/`: Key figures and visualizations
	- `abstract.png`: Visual abstract (see above)
	- `cooccurrence_heatmap.png`, `ige_expression_boxplot.png`, `LAC_cube.png`, `presence_absence_heatmap.png`, `scree_plot.png`: Analysis figures
- `requirements.txt`: Python dependencies
- `Dockerfile`: Reproducible environment setup



## How to Reproduce

1. **Clone the repository**
2. **Set up the environment**
	- Use the provided Docker image for a fully reproducible environment:
		- Pull the image: `docker pull ovlasovets/rpy2:latest`
		- Start a container (mounting your repo):
			```bash
			docker run -it --rm -p 8888:8888 -v $(pwd):/workspace ovlasovets/rpy2:latest
			```
		- Launch Jupyter inside the container as instructed, or follow the container's README for details.
	- Alternatively, install Python dependencies from `requirements.txt`:
		```bash
		pip install -r requirements.txt
		```
3. **Run the notebooks**
	- Open the notebooks in the `notebooks/` directory and follow the instructions in each notebook.


## License

MIT License


## Contact

For questions, contact oleg.vlasovets@helmholtz-munich.de.