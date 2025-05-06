# capstone
this is capstone for spring 2025

# Directions

## Python

Ensure that you have python 3.12 installed!

## Compiling the website

Run the vitual python environent in the directory capstone/
``` python
python -m venv venv
```

then on windows
``` cmd
.\venv\Scripts\activate
```

or on unix (linux, macOS)
``` Bash
source venv/bin/activate
```

**Install Dependencies**
``` python
pip install -r requirements.txt
```

*Start the app*
``` python
python run.py
```

## Working the website

Just click the run button :)


 # Deploy with Kuber

 docker build -t capstone-app:0.1 .

 kubectl apply -f deploy.yaml
kubectl apply -f svc.yaml
