FROM python:3.12.4

RUN mkdir -p /app
WORKDIR /app
ADD ./ /app/

RUN tar -xzf /app/config_vcf2circos_29032023.tar.gz -C /app/
RUN ["python", "/app/setup.py", "install"]
EXPOSE 5000

CMD ["python", "/app/main.py"]

