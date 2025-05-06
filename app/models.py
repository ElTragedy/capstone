from app import db
from datetime import datetime

class MoleculeSubmission(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    user_name = db.Column(db.String(100), nullable=False)
    score = db.Column(db.Float, nullable=False)
    smiles = db.Column(db.Text, nullable=False)
    mol_data = db.Column(db.Text, nullable=False)  
    created_at = db.Column(db.DateTime, default=datetime.utcnow)
    
    def to_dict(self):
        return {
            'id': self.id,
            'user_name': self.user_name,
            'score': self.score,
            'smiles': self.smiles,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S')
        } 