set -ev
docker_image=us.gcr.io/cds-docker-containers/depmap-crispr-vs-rnai:1
docker build -t ${docker_image} .
docker push ${docker_image}
