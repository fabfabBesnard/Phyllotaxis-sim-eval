# (Draft for a) docker readme

This readme explains:

    - how to download (pull) a pre-built docker image
    - how to run notebooks from a docker container with all the requirements (sm-dtw, R, jupyter notebooks), avoiding multiple and possibly complex install procedures
    - how to solve most current troubles
    - (advanced / for developers): how to build a new /updated docker image

## Requisite
You must install [docker](https://docs.docker.com/desktop/#download-and-install)
## 1. Download the docker image
It should be downloaded automatically when trying to run the docker for the first time
```
docker run -p 8888:8888 -it roboticsmicrofarms/sm-dtw_demo:latest bash
```

## 2. Using jupyter notebook inside a docker:
https://simplernerd.com/docker-jupyter-notebook/

### Running the Jupyter notebook:
#reference: https://www.thecodeteacher.com/question/55171/python---Access-Jupyter-notebook-running-on-Docker-container
cd path/to/gitclone/Phyllotaxis-sim-eval
docker run -p 8888:8888 -it roboticsmicrofarms/sm-dtw_demo:latest bash
#note: '/myapp' is the starting workfolder inside the Docker image

### Inside the container, serve the notebooks:
1. Enter the following command after the running container's prompt (which should be something like: `root@aa7632c1fc29:/myapp#`):
```
jupyter notebook Phyllotaxis-sim-eval/notebooks/docker_run --ip 0.0.0.0 --no-browser --allow-root
```
2. Outside the docker, start your browser
3. enter the url with token:"http://127.0.0.1:8888/?token=....." or directly Ctrl+click on it
4. select the correct kernel (R, bash or python)
5. run the notebook !

### Troobleshooting
1. Docker external connectivity
example of error message:
```
docker: Error response from daemon: driver failed programming external connectivity on endpoint vibrant_carson (91c770397f27b2a424b77a55f9253a271e64c83a1fdf898b487bfc6808193497): Error starting userland proxy: listen tcp4 0.0.0.0:8888: bind: address already in use.
ERRO[0000] error waiting for container: context canceled
```

**Cause of the problem**: another process is currently listening to the 8888 port

**solutions**: 
- Check that other running containers are not listening to the port (`docker ps`) and stop them if necessary (`` )
- find the process that is listening to the port with the command `lsof -i:8888` and consider killing that process if it is safe (`kill $(lsof -t -i:"8888")`)
- if nothing works... shut down your computer  (all process connected to the port should be stopped) and re-start !

2. *"When running the docker, R starts directly"*: you probably forgot `bash` at the end of the start command `docker run -it roboticsmicrofarms/sm-dtw_demo:latest bash`

3. Jupyter notebook connection
**Cause of the problem**: another server is running on the same port 

**solutions**: stop the other process (stopping the web-browser page is not enough, stop the process in the terminal, e.g. another jupter notebook that could be running)
## 3. Build a new image docker

### Build procedure
1. in this folder (Phyllotaxis-sim-eval/docker), just enter:
```
docker build -t imagename:tagname . #(e.g: your_docker_account/repo:latest so that you can push directly. Do not forget the '.' at the end)
```
2. Up-dating the repo `Phyllotaxis-sim-eval` in the docker build
   - for testing (e.g. developers), the simplest & fastest is to substitute your local repo of `Phyllotaxis-sim-eval` by mounting it when starting the docker container:
```
docker run -v path/to/your/local/repo:/myapp/Phyllotaxis-sim-eval -p 8888:8888 -it roboticsmicrofarms/sm-dtw_phyllotaxis:latest bash
```
### testing that everything works inside the docker
Start a docker container from the image in the interactive mode
docker run -it my_new_build:mytag bash

### References about the Dockerfile recipe used for image build:

- install R inside the docker image: https://towardsdatascience.com/creating-sandbox-environments-for-r-with-docker-def54e3491a3
the rocker project: https://github.com/rocker-org/rocker
- Installing R and bash kernels in Jupyter:
https://evodify.com/python-r-bash-jupyter-notebook/
- Install sm-dtw in the docker via miniconda
https://linuxhandbook.com/dockerize-python-apps/