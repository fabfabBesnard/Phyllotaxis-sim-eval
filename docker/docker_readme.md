# Draft for a readme

This readme explains:
    - how to build a docker image to use sm-dtw with simulated phyllotaxis data
    - how to download (pull) a pre-built docker image
    - how to run notebooks from a docker container with all the requirements (sm-dtw, R, jupyter notebooks), avoiding multiple and possibly complex install procedures

## make the Dockerfile recipe for image build:
Some reference to install R inside the docker image: https://towardsdatascience.com/creating-sandbox-environments-for-r-with-docker-def54e3491a3
the rocker project: https://github.com/rocker-org/rocker

### Installing R and bash kernels in Jupyter:
https://evodify.com/python-r-bash-jupyter-notebook/
### Install sm-dtw in teh docker via miniconda
https://linuxhandbook.com/dockerize-python-apps/


## Using jupyter notebook inside a docker:
https://simplernerd.com/docker-jupyter-notebook/

### testing that everything works inside the docker
Start a docker container from the image in the interactive mode
docker run -it physimev:test bash

### Running the Jupyter notebook:
#reference: https://www.thecodeteacher.com/question/55171/python---Access-Jupyter-notebook-running-on-Docker-container
cd path/to/gitclone/Phyllotaxis-sim-eval
docker run -v "$(pwd):/myapp/Phyllotaxis-sim-eval" -p 8888:8888 -it physimev:test bash
#note: '/myapp' is the starting workspace inside the Docker image

### Inside the container, serve the notebooks:
1. root@aa7632c1fc29:/myapp#jupyter notebook Phyllotaxis-sim-eval/notebooks --ip 0.0.0.0 --no-browser --allow-root
#Outside the docker, start your browser:
2. enter the url with token:"http://127.0.0.1:8888/?token=....." or directly click on it
3. select the correct kernel (R, bash or python)
4. run the notebook !
