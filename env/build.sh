set -ev
docker_image=jkrillbu/depmap-crispr-vs-rnai:1
docker build -t ${docker_image} .
docker push ${docker_image}
