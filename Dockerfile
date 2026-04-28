FROM python:3-slim
WORKDIR /app

# COPY requirements.txt .
# RUN pip install -r requirements.txt

COPY . .
RUN python3 -m pip install .
ENTRYPOINT ["sofia_image_pipeline"]