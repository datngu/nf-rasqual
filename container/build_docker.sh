## BUILD DOCKER AUTOMATICLY

docker build -t rasqual:v0.0.0 -f Dockerfile .

# v0.0.0
docker tag rasqual:v0.0.0 ndatth/rasqual:v0.0.0
docker push ndatth/rasqual:v0.0.0

