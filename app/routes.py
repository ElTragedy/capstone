from app import app
from flask import render_template, request, jsonify, send_file
import subprocess  # or use your own notebook runner code
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import Draw
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote, unquote
from io import BytesIO

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
    
@app.route('/get_mol_image')
def get_mol_image():
        # Read SMILES string from file
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()
    
        # Compute the 2D coordinates of the molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({"error": "Invalid SMILES string"}), 400
        
        AllChem.Compute2DCoords(mol)
        
        # Generate the image
        img = Draw.MolToImage(mol)
        
        # Save the image to a byte stream
        img_byte_array = BytesIO()
        img.save(img_byte_array, format='PNG')
        img_byte_array.seek(0)
        
        # Return the image as a response
        return send_file(img_byte_array, mimetype='image/png')
        

@app.route('/get_mol_info')
def get_mol_info():
    try:
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()

        print(smiles)
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' + quote(smiles) + '/property/Title,MolecularFormula/JSON'
        print(url)
        response = requests.get(url)
        data = response.json()
        print(data)
            
        return data
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/get_mol_creation_date')
def get_mol_creation_date():
    try:
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()

        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/' + quote(smiles) + '/dates/JSON'
        response = requests.get(url)
        data = response.json()
            
        return data
    except Exception as e:
        return jsonify({"error": str(e)}), 500