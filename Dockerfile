FROM python:3.8-slim
WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY src/ .
ENTRYPOINT ["python3", "image_pipeline.py"]