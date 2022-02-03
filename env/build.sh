set -ev
docker_image=us.gcr.io/cds-docker-containers/depmap-crispr-vs-rnai:1
docker build -f env/Dockerfile -t ${docker_image} .
docker push ${docker_image}
