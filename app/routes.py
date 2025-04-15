from app import app, db
from flask import render_template, request, jsonify, send_file, redirect, url_for
import subprocess  # or use your own notebook runner code
import sys
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote, unquote
from io import BytesIO
from app.models import MoleculeSubmission
import random  # For random score generation


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/run_notebook", methods=["POST"])
def run_notebook():
    try:
        # Run nbconvert to generate output.html from the notebook
        subprocess.check_call(
            [
                sys.executable,
                "-m",
                "nbconvert",
                "--to",
                "html",
                "--execute",
                "notebooks/smiles_drawing.ipynb",
                "--output",
                "output.html",
                "--output-dir",
                ".",
            ]
        )
        with open("output.html", "r") as f:
            output = f.read()

        # Use BeautifulSoup to extract only the body content
        soup = BeautifulSoup(output, "html.parser")
        body_content = soup.body.decode_contents()

        # Generate a random score between 0 and 100
        # TODO: Replace this with actual scoring logic based on molecule properties
        score = round(random.uniform(60, 100), 2)

        return jsonify({"output": body_content, "score": score})
    except subprocess.CalledProcessError as e:
        return jsonify({"error": str(e)}), 500


@app.route("/publish_molecule", methods=["POST"])
def publish_molecule():
    try:
        data = request.json
        user_name = data.get('user_name')
        score = data.get('score')
        
        # Read the current molecule data
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()
        
        # Get MOL data
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mol_block = Chem.MolToMolBlock(mol)
        
        # Create new submission
        submission = MoleculeSubmission(
            user_name=user_name,
            score=score,
            smiles=smiles,
            mol_data=mol_block
        )
        
        db.session.add(submission)
        db.session.commit()
        
        return jsonify({"success": True, "submission_id": submission.id})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/leaderboard")
def leaderboard():
    # Get top submissions, ordered by score
    submissions = MoleculeSubmission.query.order_by(MoleculeSubmission.score.desc()).all()
    return render_template("leaderboard.html", submissions=submissions)


@app.route("/view_molecule/<int:submission_id>")
def view_molecule(submission_id):
    submission = MoleculeSubmission.query.get_or_404(submission_id)
    
    # Write the SMILES to the file so the existing viewer can use it
    with open("notebooks/molecule.smiles", "w") as file:
        file.write(submission.smiles)
    
    return render_template("view_molecule.html", submission=submission)


@app.route("/get_molecule")
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


@app.route("/get_mol_image")
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
    width, height = 300, 300
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)

    # Draw options
    drawer.drawOptions().clearBackground = (
        False  # This makes the background transparent
    )

    # Draw the molecule
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    # Get PNG image bytes
    png_data = drawer.GetDrawingText()

    # Return as a response
    return send_file(BytesIO(png_data), mimetype="image/png")


@app.route("/get_mol_info")
def get_mol_info():
    try:
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()

        print(smiles)
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            + quote(smiles)
            + "/property/Title,MolecularFormula/JSON"
        )
        print(url)
        response = requests.get(url)
        data = response.json()
        
        # Add style to the title to make it white
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable'] and len(data['PropertyTable']['Properties']) > 0:
            title = data['PropertyTable']['Properties'][0].get('Title', '')
            data['PropertyTable']['Properties'][0]['Title'] = f'<span style="color: #ecf0f1;">{title}</span>'
        
        print(data)
        return data
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/get_mol_creation_date")
def get_mol_creation_date():
    try:
        with open("notebooks/molecule.smiles", "r") as file:
            smiles = file.read().strip()

        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
            + quote(smiles)
            + "/dates/JSON"
        )
        response = requests.get(url)
        data = response.json()

        return data
    except Exception as e:
        return jsonify({"error": str(e)}), 500
