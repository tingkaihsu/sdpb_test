# Usage of Mathematica code
First one must have mathematica or wolfram engine, and note that the provide files merely define the function that can be used.
One can modify the `.m` content or call the function by:
```bash
wolframscript -file filename
```
# Using `Docker` to run `SDPB`
Suppose you have an input file `/my/project/input.json`. To use this
file, and be able to write output, we will make the `/my/project`
directory visible to the docker image in the location `/usr/local/share/sdpb`.

The structure of the [docker run](https://docs.docker.com/engine/reference/commandline/run/) command is

    docker run <options> <image> <command>
    
In this case, `<options>` will be used to mount your directory in a
place that docker will see it.

    -v /my/project/:/usr/local/share/sdpb/

`<image>` is the image name, e.g. `bootstrapcollaboration/sdpb:master`

`<command>` is the command that you would normally use to run the SDPB
commands (see [Usage.md](Usage.md)).  The directory containing the
input file is mounted as `/usr/local/share/sdpb`. So we first run `pmp2sdp` to
convert from json

    mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/input.json -o /usr/local/share/sdpb/sdp
    
`mpirun` runs as root inside the docker container.  Running `mpirun` as
root is normally dangerous, but it is safe to do so inside the
container. To allow `mpirun` to run as root, we add the option
`--allow-run-as-root`. This uses 4 cores when running pmp2sdp. You
can change that number to match your own machine.  Putting it all
together on a single line

    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 pmp2sdp --precision 1024 -i /usr/local/share/sdpb/input.json -o /usr/local/share/sdpb/sdp

Running this command will create directory `/my/project/sdp`.
To search for primal-dual solutions

    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master mpirun --allow-run-as-root -n 4 sdpb --precision=1024 -s /usr/local/share/sdpb/sdp

The results will be in `/my/project/`.

Note that the newly created files may be owned by root.
If you cannot remove them outside the container, run `rm` from the container, e.g.:

    docker run -v /my/project/:/usr/local/share/sdpb bootstrapcollaboration/sdpb:master rm /usr/local/share/sdpb/sdp
