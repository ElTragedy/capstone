FROM tensorflow/tensorflow:latest-gpu

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# ✅ Add PyFiles to the PYTHONPATH
#ENV PYTHONPATH="${PYTHONPATH}:/app/PyFiles"

WORKDIR /app

COPY requirements.txt /app/
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY . /app/


CMD ["python", "train.py"]
