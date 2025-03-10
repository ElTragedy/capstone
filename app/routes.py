from app import app
from flask import render_template, request, jsonify
import subprocess  # or use your own notebook runner code
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import requests
from bs4 import BeautifulSoup

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/run_notebook', methods=['POST'])
def run_notebook():
    try:
        # Example using nbconvert to run and convert your notebook to HTML:
        subprocess.check_call([
            sys.executable, '-m', 'nbconvert',
            '--to', 'html',
            '--execute', 'notebooks/smiles_drawing.ipynb',
            '--output', 'output.html',
            '--output-dir', '.'
        ])
        with open('output.html', 'r') as f:
            output = f.read()
        return jsonify({'output': output})
    except subprocess.CalledProcessError as e:
        return jsonify({'error': str(e)}), 500

@app.route('/get_molecule')
def get_molecule():
    try:
        # Read SMILES string from file
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()

        # Convert SMILES to 3D molecule
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)  # Add hydrogen atoms
        AllChem.EmbedMolecule(mol)  # Generate 3D coordinates

        # Convert to MOL format
        mol_block = Chem.MolToMolBlock(mol)
        
        print(mol_block)
        
        return mol_block
    except Exception as e:
        return jsonify({"error": str(e)}), 500

    @app.route('/get_mol_info')
    def get_mol_info():
        try:
            #read from this link https://pubchem.ncbi.nlm.nih.gov/#query=[smile here]
            #scrape the data and return it
            
            smile = 'notebooks/molecule.smiles'
            url = 'https://pubchem.ncbi.nlm.nih.gov/#query=' + smile
            page = requests.get(url)
            soup = BeautifulSoup(page.content, 'html.parser')
            mol_info = soup.find_all('div', class_='f-0875')
            print(mol_info)
            
            return jsonify(mol_info)
        except Exception as e:
            return jsonify({"error": str(e)}), 500