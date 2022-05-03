FROM python:3.8-slim
WORKDIR /app

COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .
RUN python3 setup.py install
ENTRYPOINT ["sofia_image_pipeline"]