# Slim Python 3.11 image
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY requirements.txt .

# Install any needed packages
RUN apt-get update && apt-get install -y \
    gcc \
    clustalo 

RUN pip3 install --no-cache-dir numpy biopython

# Copy the rest of your application's source code into the container
COPY . .

# Run the application
CMD ["python", "main.py"]
