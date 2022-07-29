FROM python:latest

WORKDIR /app

# By copying over requirements first, we make sure that Docker will cache
# our installed requirements rather than reinstall them on every build
COPY requirements.txt /app/
RUN pip install --no-cache-dir -r /app/requirements.txt