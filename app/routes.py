from app import app
from flask import render_template, request, jsonify
import subprocess  # or use your own notebook runner code
import sys

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
